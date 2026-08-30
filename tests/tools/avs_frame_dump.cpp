#include <avisynth.h>

#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <dlfcn.h>
#endif

const AVS_Linkage* AVS_linkage = nullptr;

namespace {

constexpr int kMinimumAviSynthInterfaceVersion = 6;

struct Options {
  std::filesystem::path script;
  std::filesystem::path raw_out;
  std::filesystem::path meta_out;
  std::string runtime;
  std::vector<std::pair<std::string, int>> expected_int_vars;
  int frame = 0;
};

class DynamicLibrary {
public:
  DynamicLibrary() = default;
  DynamicLibrary(const DynamicLibrary&) = delete;
  DynamicLibrary& operator=(const DynamicLibrary&) = delete;

  DynamicLibrary(DynamicLibrary&& other) noexcept
    : handle_(other.handle_) {
    other.handle_ = nullptr;
  }

  DynamicLibrary& operator=(DynamicLibrary&& other) noexcept {
    if (this != &other) {
      close();
      handle_ = other.handle_;
      other.handle_ = nullptr;
    }
    return *this;
  }

  ~DynamicLibrary() {
    close();
  }

  static DynamicLibrary open(const std::string& path) {
    DynamicLibrary library;
#if defined(_WIN32)
    library.handle_ = LoadLibraryA(path.c_str());
    if (library.handle_ == nullptr) {
      throw std::runtime_error("failed to load " + path);
    }
#else
    library.handle_ = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
    if (library.handle_ == nullptr) {
      const char* error = dlerror();
      throw std::runtime_error(
        "failed to load " + path + ": " + (error != nullptr ? error : "unknown dlopen error")
      );
    }
#endif
    return library;
  }

  template <class Function>
  Function symbol(const char* name) const {
#if defined(_WIN32)
    auto* pointer = reinterpret_cast<Function>(GetProcAddress(static_cast<HMODULE>(handle_), name));
#else
    dlerror();
    auto* pointer = reinterpret_cast<Function>(dlsym(handle_, name));
#endif
    if (pointer == nullptr) {
#if defined(_WIN32)
      throw std::runtime_error(std::string{"failed to load symbol "} + name);
#else
      const char* error = dlerror();
      throw std::runtime_error(
        std::string{"failed to load symbol "} + name + ": " +
        (error != nullptr ? error : "unknown dlsym error")
      );
#endif
    }
    return pointer;
  }

private:
  void close() {
    if (handle_ == nullptr) {
      return;
    }
#if defined(_WIN32)
    FreeLibrary(static_cast<HMODULE>(handle_));
#else
    dlclose(handle_);
#endif
    handle_ = nullptr;
  }

#if defined(_WIN32)
  HMODULE handle_ = nullptr;
#else
  void* handle_ = nullptr;
#endif
};

using CreateScriptEnvironmentFn = IScriptEnvironment*(__stdcall *)(int);

[[noreturn]] void throw_usage() {
  throw std::runtime_error(
    "usage: f3kdb_avs_frame_dump --script script.avs --frame n "
    "--raw-out frame.bin --meta-out frame.json [--runtime path] "
    "[--expect-int-var name=value ...]"
  );
}

int parse_int(std::string_view text, const char* name) {
  std::size_t used = 0;
  int value = 0;
  try {
    value = std::stoi(std::string{text}, &used, 10);
  } catch (const std::exception&) {
    throw std::runtime_error(std::string{"invalid "} + name + ": " + std::string{text});
  }
  if (used != text.size()) {
    throw std::runtime_error(std::string{"invalid "} + name + ": " + std::string{text});
  }
  return value;
}

std::pair<std::string, int> parse_expected_int_var(std::string_view text) {
  const std::size_t separator = text.find('=');
  if (separator == std::string_view::npos || separator == 0 || separator + 1 == text.size()) {
    throw std::runtime_error("invalid expected int variable: " + std::string{text});
  }

  std::string name{text.substr(0, separator)};
  const int value = parse_int(text.substr(separator + 1), "expected int variable value");
  return {std::move(name), value};
}

std::string_view require_next_arg(int& index, int argc, char** argv) {
  if (++index >= argc) {
    throw_usage();
  }
  return argv[index];
}

Options parse_args(int argc, char** argv) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string_view arg = argv[i];
    if (arg == "--script") {
      options.script = require_next_arg(i, argc, argv);
    } else if (arg == "--frame") {
      options.frame = parse_int(require_next_arg(i, argc, argv), "frame");
    } else if (arg == "--raw-out") {
      options.raw_out = require_next_arg(i, argc, argv);
    } else if (arg == "--meta-out") {
      options.meta_out = require_next_arg(i, argc, argv);
    } else if (arg == "--runtime") {
      options.runtime = require_next_arg(i, argc, argv);
    } else if (arg == "--expect-int-var") {
      options.expected_int_vars.push_back(parse_expected_int_var(require_next_arg(i, argc, argv)));
    } else {
      throw_usage();
    }
  }

  if (options.script.empty() || options.raw_out.empty() || options.meta_out.empty()) {
    throw_usage();
  }
  if (options.frame < 0) {
    throw std::runtime_error("frame must be non-negative");
  }
  return options;
}

std::vector<std::string> runtime_candidates(const Options& options) {
  if (!options.runtime.empty()) {
    return {options.runtime};
  }

#if defined(_WIN32)
  return {"avisynth.dll"};
#elif defined(__APPLE__)
  return {"libavisynth.dylib", "libavisynth.so", "libavisynth.so.11"};
#else
  return {"libavisynth.so", "libavisynth.so.11"};
#endif
}

DynamicLibrary load_avisynth_runtime(const Options& options) {
  std::string errors;
  for (const auto& candidate : runtime_candidates(options)) {
    try {
      return DynamicLibrary::open(candidate);
    } catch (const std::exception& error) {
      if (!errors.empty()) {
        errors += "; ";
      }
      errors += error.what();
    }
  }
  throw std::runtime_error("could not load AviSynth+ runtime: " + errors);
}

IScriptEnvironment* create_script_environment(CreateScriptEnvironmentFn create_environment) {
  std::string attempted_versions;
  for (
    int version = AVISYNTH_INTERFACE_VERSION;
    version >= kMinimumAviSynthInterfaceVersion;
    --version
  ) {
    if (!attempted_versions.empty()) {
      attempted_versions += ", ";
    }
    attempted_versions += std::to_string(version);
    if (IScriptEnvironment* env = create_environment(version)) {
      return env;
    }
  }

  throw std::runtime_error(
    "CreateScriptEnvironment returned null for supported interface versions: " + attempted_versions
  );
}

PClip import_clip(IScriptEnvironment* env, const std::filesystem::path& script) {
  const auto absolute_script = std::filesystem::absolute(script);
  const auto parent = absolute_script.parent_path();
  if (!parent.empty()) {
    const auto parent_string = parent.string();
    env->SetWorkingDir(parent_string.c_str());
  }

  const auto script_string = absolute_script.string();
  AVSValue arg(script_string.c_str());
  AVSValue result = env->Invoke("Import", AVSValue(&arg, 1));
  if (!result.IsClip()) {
    throw std::runtime_error("script did not return an AviSynth clip: " + script_string);
  }
  return result.AsClip();
}

void verify_expected_int_vars(
  IScriptEnvironment* env,
  const std::vector<std::pair<std::string, int>>& expected_int_vars
) {
  for (const auto& [name, expected] : expected_int_vars) {
    AVSValue value = env->GetVarDef(name.c_str(), AVSValue());
    if (!value.Defined()) {
      throw std::runtime_error("AviSynth variable was not set: " + name);
    }
    if (!value.IsInt()) {
      throw std::runtime_error("AviSynth variable is not an int: " + name);
    }
    const int actual = value.AsInt();
    if (actual != expected) {
      throw std::runtime_error(
        "AviSynth variable has unexpected value: " + name +
        " expected " + std::to_string(expected) +
        " got " + std::to_string(actual)
      );
    }
  }
}

std::vector<std::pair<int, const char*>> canonical_planes(const VideoInfo& vi) {
  if (vi.IsPlanarRGBA()) {
    return {{PLANAR_R, "R"}, {PLANAR_G, "G"}, {PLANAR_B, "B"}, {PLANAR_A, "A"}};
  }
  if (vi.IsPlanarRGB()) {
    return {{PLANAR_R, "R"}, {PLANAR_G, "G"}, {PLANAR_B, "B"}};
  }
  if (vi.IsYUVA()) {
    return {{PLANAR_Y, "Y"}, {PLANAR_U, "U"}, {PLANAR_V, "V"}, {PLANAR_A, "A"}};
  }
  if (vi.IsY()) {
    return {{PLANAR_Y, "Y"}};
  }
  if (vi.IsYUV()) {
    return {{PLANAR_Y, "Y"}, {PLANAR_U, "U"}, {PLANAR_V, "V"}};
  }

  throw std::runtime_error("unsupported AviSynth output format for baseline hashing");
}

void write_frame_data(
  const Options& options,
  const VideoInfo& vi,
  const PVideoFrame& frame
) {
  std::ofstream raw(options.raw_out, std::ios::binary);
  if (!raw) {
    throw std::runtime_error("failed to open raw output: " + options.raw_out.string());
  }

  const auto planes = canonical_planes(vi);
  for (const auto& [plane, name] : planes) {
    (void)name;
    const BYTE* src = frame->GetReadPtr(plane);
    const int pitch = frame->GetPitch(plane);
    const int row_size = frame->GetRowSize(plane);
    const int height = frame->GetHeight(plane);
    for (int y = 0; y < height; ++y) {
      const BYTE* row = src + (static_cast<std::ptrdiff_t>(y) * pitch);
      raw.write(reinterpret_cast<const char*>(row), row_size);
    }
  }

  std::ofstream meta(options.meta_out);
  if (!meta) {
    throw std::runtime_error("failed to open metadata output: " + options.meta_out.string());
  }

  meta << "{\n";
  meta << "  \"width\": " << vi.width << ",\n";
  meta << "  \"height\": " << vi.height << ",\n";
  meta << "  \"frames\": " << vi.num_frames << ",\n";
  meta << "  \"pixel_type\": " << vi.pixel_type << ",\n";
  meta << "  \"bits_per_sample\": " << vi.BitsPerComponent() << ",\n";
  meta << "  \"bytes_per_sample\": " << vi.ComponentSize() << ",\n";
  meta << "  \"planes\": [\n";
  for (std::size_t i = 0; i < planes.size(); ++i) {
    const auto& [plane, name] = planes[i];
    meta << "    {"
         << R"("name": ")" << name << "\", "
         << "\"row_bytes\": " << frame->GetRowSize(plane) << ", "
         << "\"height\": " << frame->GetHeight(plane) << ", "
         << "\"stride\": " << frame->GetPitch(plane) << "}";
    if (i + 1 != planes.size()) {
      meta << ",";
    }
    meta << "\n";
  }
  meta << "  ]\n";
  meta << "}\n";
}

void run(const Options& options, IScriptEnvironment* env) {
  PClip clip = import_clip(env, options.script);
  verify_expected_int_vars(env, options.expected_int_vars);
  const VideoInfo& vi = clip->GetVideoInfo();
  if (!vi.HasVideo()) {
    throw std::runtime_error("script returned a clip without video");
  }
  if (options.frame >= vi.num_frames) {
    throw std::runtime_error("requested frame is outside the clip frame range");
  }

  PVideoFrame frame = clip->GetFrame(options.frame, env);
  if (!frame) {
    throw std::runtime_error("GetFrame returned an empty frame");
  }
  write_frame_data(options, vi, frame);
}

void run_with_environment(const Options& options, CreateScriptEnvironmentFn create_environment) {
  IScriptEnvironment* env = create_script_environment(create_environment);

  try {
    AVS_linkage = env->GetAVSLinkage();
    run(options, env);
  } catch (const AvisynthError& error) {
    const std::string message = error.msg != nullptr ? error.msg : "unknown AviSynth error";
    env->DeleteScriptEnvironment();
    AVS_linkage = nullptr;
    throw std::runtime_error("AviSynth error: " + message);
  } catch (...) {
    env->DeleteScriptEnvironment();
    AVS_linkage = nullptr;
    throw;
  }

  env->DeleteScriptEnvironment();
  AVS_linkage = nullptr;
}

} // namespace

int main(int argc, char** argv) {
  try {
    const Options options = parse_args(argc, argv);
    DynamicLibrary runtime = load_avisynth_runtime(options);
    const auto create_environment =
      runtime.symbol<CreateScriptEnvironmentFn>("CreateScriptEnvironment");
    run_with_environment(options, create_environment);
    return EXIT_SUCCESS;
  } catch (const AvisynthError& error) {
    std::cerr << "AviSynth error: " << error.msg << '\n';
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
  }
  return EXIT_FAILURE;
}
