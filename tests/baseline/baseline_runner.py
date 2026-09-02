from __future__ import annotations

import argparse
import ctypes
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Any

GOLDEN_SCHEMA = 2
PLUGIN_NAMESPACE = "neo_f3kdb"
DEFAULT_FILTER = "Deband"
BASELINE_DIR = Path(__file__).resolve().parent

class PlaneBytes:
    def __init__(self, name: str, width: int, height: int, bytes_per_sample: int, stride: int, data: bytes):
        self.name = name
        self.width = width
        self.height = height
        self.bytes_per_sample = bytes_per_sample
        self.stride = stride
        self.data = data

    @property
    def row_bytes(self) -> int:
        return self.width * self.bytes_per_sample

def hash_planes(planes: Iterable[PlaneBytes]) -> str:
    digest = hashlib.sha256()
    for plane in planes:
        row_bytes = plane.row_bytes
        for y in range(plane.height):
            start = y * plane.stride
            digest.update(plane.data[start:start + row_bytes])
    return digest.hexdigest()

def select_cases(cases: list[dict], tier: str) -> list[dict]:
    if tier == "all":
        return list(cases)
    return [case for case in cases if case.get("tier") == tier]

def load_cases(cases_dir: Path) -> list[dict]:
    cases: list[dict] = []
    for path in sorted(cases_dir.glob("*.json")):
        with path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        cases.extend(_resolve_case_sources(payload.get("cases", []), path.parent))
    return cases

def _resolve_case_sources(cases: list[dict], base_dir: Path) -> list[dict]:
    resolved_cases = []
    for case in cases:
        resolved_case = dict(case)
        source = dict(case["source"])
        if source["type"] in ("ffms2", "bestsource"):
            source["resolved_path"] = (base_dir / source["path"]).resolve()
        resolved_case["source"] = source
        resolved_cases.append(resolved_case)
    return resolved_cases

_LOADED_VS_PLUGINS: set[str] = set()

def run_vs_case(case: dict, plugin_path: Path, backend: str) -> list[dict]:
    import vapoursynth as vs

    if os.name == "nt":
        os.add_dll_directory(str(plugin_path.parent))

    core = vs.core
    plugin_key = str(plugin_path.resolve())
    if plugin_key not in _LOADED_VS_PLUGINS:
        core.std.LoadPlugin(path=plugin_key)
        _LOADED_VS_PLUGINS.add(plugin_key)

    source = case["source"]
    clip = _vs_source_clip(core, vs, source)
    
    params = dict(case.get("params", {}))
    if backend == "purec":
        params["opt"] = 0
    elif backend == "highway":
        params.pop("opt", None)

    filter_name = case.get("filter", DEFAULT_FILTER)
    namespace = getattr(core, PLUGIN_NAMESPACE)
    clip = getattr(namespace, filter_name)(clip, **params)

    results = []
    for frame_number in case["frames"]:
        frame = clip.get_frame(frame_number)
        planes, metadata = _vs_frame_planes(frame, clip.num_frames)
        results.append(_result_record(case, "vs", frame_number, hash_planes(planes), metadata, backend))
    return results

def _vs_source_clip(core: Any, vs: Any, source: dict) -> Any:
    if source["type"] == "blank":
        return core.std.BlankClip(
            width=source["width"],
            height=source["height"],
            length=source["length"],
            format=getattr(vs, source["format"]),
            color=source["color"],
        )
    if source["type"] == "blank_sequence":
        clips = [
            core.std.BlankClip(
                width=source["width"],
                height=source["height"],
                length=1,
                format=getattr(vs, source["format"]),
                color=color,
            )
            for color in source["colors"]
        ]
        return core.std.Splice(clips)
    if source["type"] == "ffms2":
        path = str(source["resolved_path"])
        if hasattr(core, "ffms2"):
            return core.ffms2.Source(source=path)
        if hasattr(core, "lsmas"):
            return core.lsmas.LibavSMASHSource(source=path)
        raise RuntimeError(
            "no supported VapourSynth source filter found; tried "
            "core.ffms2.Source and core.lsmas.LibavSMASHSource"
        )
    if source["type"] == "synthetic_gradient":
        blank = core.std.BlankClip(
            width=source["width"],
            height=source["height"],
            length=source["length"],
            format=getattr(vs, source["format"])
        )
        def populate_gradient(n, f):
            import numpy as np
            f_out = f.copy()
            for plane_index in range(f_out.format.num_planes):
                sub_w = f_out.format.subsampling_w if plane_index > 0 else 0
                sub_h = f_out.format.subsampling_h if plane_index > 0 else 0
                w = f_out.width >> sub_w
                h = f_out.height >> sub_h
                
                y_idx, x_idx = np.ogrid[:h, :w]
                
                if plane_index == 0:
                    val = x_idx + y_idx
                    if "10" in source["format"]:
                        val = val * 4
                    elif "12" in source["format"]:
                        val = val * 16
                    elif "16" in source["format"]:
                        val = val * 256
                elif plane_index == 1:
                    if "10" in source["format"]:
                        val = (x_idx - y_idx) * 4 + 512
                    elif "12" in source["format"]:
                        val = (x_idx - y_idx) * 16 + 2048
                    elif "16" in source["format"]:
                        val = (x_idx - y_idx) * 256 + 32768
                    else:
                        val = x_idx - y_idx + 128
                elif plane_index == 2:
                    if "10" in source["format"]:
                        val = (y_idx - x_idx) * 4 + 512
                    elif "12" in source["format"]:
                        val = (y_idx - x_idx) * 16 + 2048
                    elif "16" in source["format"]:
                        val = (y_idx - x_idx) * 256 + 32768
                    else:
                        val = y_idx - x_idx + 128
                        
                bits = f_out.format.bits_per_sample
                max_val = (1 << bits) - 1
                val = np.clip(val, 0, max_val)
                
                p = f_out.get_write_ptr(plane_index)
                stride = f_out.get_stride(plane_index)
                
                dtype = np.uint16 if f_out.format.bytes_per_sample == 2 else np.uint8
                ptr_val = ctypes.cast(p, ctypes.c_void_p).value
                arr = np.ctypeslib.as_array(ctypes.cast(ptr_val, ctypes.POINTER(ctypes.c_uint16 if dtype == np.uint16 else ctypes.c_uint8)), shape=(h, stride // f_out.format.bytes_per_sample))
                arr[:, :w] = val
                
            return f_out
            
        return core.std.ModifyFrame(blank, blank, populate_gradient)
    if source["type"] == "bestsource":
        path = str(source["resolved_path"])
        clip = None
        if hasattr(core, "bs"):
            clip = core.bs.VideoSource(source=path)
        elif hasattr(core, "ffms2"):
            clip = core.ffms2.Source(source=path)
        elif hasattr(core, "lsmas"):
            clip = core.lsmas.LibavSMASHSource(source=path)
        else:
            raise RuntimeError(
                "bestsource is requested but not found, and no supported VapourSynth fallback source filter was found; "
                "tried core.ffms2.Source and core.lsmas.LibavSMASHSource"
            )
        if source.get("gray"):
            clip = core.std.ShufflePlanes(clips=clip, planes=0, colorfamily=vs.GRAY)
        return clip
    raise ValueError(f"unsupported source type: {source['type']}")

def _vs_frame_planes(frame: Any, num_frames: int) -> tuple[list[PlaneBytes], dict]:
    import vapoursynth as vs
    fmt = frame.format
    planes: list[PlaneBytes] = []
    plane_metadata = []
    names = ["Y", "U", "V"] if fmt.color_family == vs.YUV else ["Y"]

    for plane_index, name in enumerate(names[:fmt.num_planes]):
        sub_w = fmt.subsampling_w if plane_index > 0 else 0
        sub_h = fmt.subsampling_h if plane_index > 0 else 0
        width = frame.width >> sub_w
        height = frame.height >> sub_h
        stride = frame.get_stride(plane_index)
        row_bytes = width * fmt.bytes_per_sample
        pointer = frame.get_read_ptr(plane_index)
        ptr_val = int(pointer) if hasattr(pointer, "__int__") else (pointer.value if hasattr(pointer, "value") else pointer)
        data = ctypes.string_at(ptr_val, stride * height)
        planes.append(PlaneBytes(name, width, height, fmt.bytes_per_sample, stride, data))
        plane_metadata.append(
            {
                "name": name,
                "row_bytes": row_bytes,
                "height": height,
                "stride": stride,
            }
        )

    metadata = {
        "width": frame.width,
        "height": frame.height,
        "frames": num_frames,
        "format": {
            "name": fmt.name,
            "id": fmt.id,
            "bits_per_sample": fmt.bits_per_sample,
            "bytes_per_sample": fmt.bytes_per_sample,
            "color_family": int(fmt.color_family),
            "sample_type": int(fmt.sample_type),
            "num_planes": fmt.num_planes,
            "subsampling_w": fmt.subsampling_w,
            "subsampling_h": fmt.subsampling_h,
        },
        "planes": plane_metadata,
    }
    return planes, metadata

def _render_avs_source_lines(source: dict) -> list[str]:
    if source["type"] == "blank_sequence":
        pixel_type = source.get("avs_format", source["format"])
        if pixel_type == "YUV420P8":
            pixel_type = "YV12"
        clips = []
        for color in source["colors"]:
            y, u, v = color
            if "10" in source["format"]:
                y >>= 2; u >>= 2; v >>= 2
            elif "12" in source["format"]:
                y >>= 4; u >>= 4; v >>= 4
            elif "16" in source["format"]:
                y >>= 8; u >>= 8; v >>= 8
            clips.append(
                f'BlankClip(width={source["width"]}, height={source["height"]}, '
                f'length=1, pixel_type="{pixel_type}", color_yuv=${y:02X}{u:02X}{v:02X})'
            )
        return ["src = " + " ++ ".join(clips)]
    elif source["type"] == "blank":
        pixel_type = source.get("avs_format", source["format"])
        if pixel_type == "YUV420P8":
            pixel_type = "YV12"
        y, u, v = source["color"]
        if "10" in source["format"]:
            y >>= 2; u >>= 2; v >>= 2
        elif "12" in source["format"]:
            y >>= 4; u >>= 4; v >>= 4
        elif "16" in source["format"]:
            y >>= 8; u >>= 8; v >>= 8
        return [
            f'src = BlankClip(width={source["width"]}, height={source["height"]}, '
            f'length={source["length"]}, pixel_type="{pixel_type}", color_yuv=${y:02X}{u:02X}{v:02X})'
        ]
    elif source["type"] == "ffms2":
        path = str(source["resolved_path"]).replace('\\', '/')
        return [
            'src = FunctionExists("FFVideoSource") '
            f'? FFVideoSource("{path}") '
            f': FunctionExists("LSMASHVideoSource") '
            f'? LSMASHVideoSource("{path}") '
            ': Assert(false, "no supported AviSynth source filter found; '
            'tried FFVideoSource and LSMASHVideoSource")'
        ]
    elif source["type"] == "synthetic_gradient":
        pixel_type = source.get("avs_format", source["format"])
        if pixel_type == "YUV420P8":
            pixel_type = "YV12"
        elif pixel_type == "GRAY8":
            pixel_type = "Y8"
            
        expr_y = "sx sy +"
        
        if "10" in source["format"]:
            expr_y = f"{expr_y} 4 *"
            expr_u = "sx sy - 4 * 512 +"
            expr_v = "sy sx - 4 * 512 +"
        elif "12" in source["format"]:
            expr_y = f"{expr_y} 16 *"
            expr_u = "sx sy - 16 * 2048 +"
            expr_v = "sy sx - 16 * 2048 +"
        elif "16" in source["format"]:
            expr_y = f"{expr_y} 256 *"
            expr_u = "sx sy - 256 * 32768 +"
            expr_v = "sy sx - 256 * 32768 +"
        else:
            expr_u = "sx sy - 128 +"
            expr_v = "sy sx - 128 +"
            
        if "GRAY" in source["format"] or pixel_type in ("Y8", "Y", "Y10", "Y12", "Y14", "Y16"):
            return [
                f'blank = BlankClip(width={source["width"]}, height={source["height"]}, '
                f'length={source["length"]}, pixel_type="{pixel_type}", color_yuv=$808080)',
                f'src = Expr(blank, "{expr_y}")'
            ]
        return [
            f'blank = BlankClip(width={source["width"]}, height={source["height"]}, '
            f'length={source["length"]}, pixel_type="{pixel_type}", color_yuv=$808080)',
            f'src = Expr(blank, "{expr_y}", "{expr_u}", "{expr_v}")'
        ]
    elif source["type"] == "bestsource":
        path = str(source["resolved_path"]).replace('\\', '/')
        lines = [
            'src = FunctionExists("BSVideoSource") '
            f'? BSVideoSource("{path}") '
            ': FunctionExists("FFVideoSource") '
            f'? FFVideoSource("{path}") '
            ': FunctionExists("LSMASHVideoSource") '
            f'? LSMASHVideoSource("{path}") '
            ': FunctionExists("LWLibavVideoSource") '
            f'? LWLibavVideoSource("{path}") '
            ': Assert(false, "no supported AviSynth source filter found; '
            'tried BSVideoSource, FFVideoSource, LSMASHVideoSource and LWLibavVideoSource")'
        ]
        if source.get("gray"):
            lines.append("src = ConvertToY(src)")
        return lines
    raise ValueError(f"unsupported source type: {source['type']}")

def run_avs_case(case: dict, plugin_path: Path, avs_dump: Path, backend: str) -> list[dict]:
    source = case["source"]
    params = dict(case.get("params", {}))
    if backend == "purec":
        params["opt"] = 0
    elif backend == "highway":
        params.pop("opt", None)

    filter_name = case.get("filter", DEFAULT_FILTER)
    avs_filter_name = "neo_f3kdb" if filter_name == "Deband" else filter_name

    param_items = []
    for key, val in sorted(params.items()):
        if isinstance(val, bool):
            val_str = "true" if val else "false"
        elif isinstance(val, str):
            val_str = f'"{val}"'
        else:
            val_str = str(val)
        param_items.append(f"{key}={val_str}")
        
    param_text = ", ".join(param_items)
    if param_text:
        param_text = ", " + param_text

    source_lines = _render_avs_source_lines(source)
    script = "\n".join(
        [
            f'LoadPlugin("{plugin_path.resolve().as_posix()}")',
            *source_lines,
            f'return {avs_filter_name}(src{param_text})',
            "",
        ]
    )

    results = []
    with tempfile.TemporaryDirectory(prefix="neo-f3kdb-baseline-") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        script_path = temp_dir / f"{case['id']}.avs"
        script_path.write_text(script, encoding="utf-8")

        for frame_number in case["frames"]:
            raw_path = temp_dir / f"{case['id']}-{frame_number}.bin"
            meta_path = temp_dir / f"{case['id']}-{frame_number}.json"
            subprocess.run(
                [
                    str(avs_dump),
                    "--script",
                    str(script_path),
                    "--frame",
                    str(frame_number),
                    "--raw-out",
                    str(raw_path),
                    "--meta-out",
                    str(meta_path),
                ],
                check=True,
            )
            sha256 = hashlib.sha256(raw_path.read_bytes()).hexdigest()
            metadata = json.loads(meta_path.read_text(encoding="utf-8"))
            results.append(_result_record(case, "avs", frame_number, sha256, metadata, backend))
    return results

def _result_record(case: dict, host: str, frame_number: int, sha256: str, metadata: dict, backend: str) -> dict:
    return {
        "case": case["id"],
        "tier": case["tier"],
        "host": host,
        "frame": frame_number,
        "source": case["source"]["id"],
        "backend": backend,
        "params": dict(case.get("params", {})),
        "sha256": sha256,
        "metadata": metadata,
    }

def write_golden(path: Path, results: list[dict]) -> None:
    payload = {"schema": GOLDEN_SCHEMA, "results": results}
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")

def verify_golden(path: Path, results: list[dict], hosts: list[str] | None = None) -> None:
    expected = json.loads(path.read_text(encoding="utf-8"))
    actual = {"schema": GOLDEN_SCHEMA, "results": results}
    
    expected_results = expected["results"]
    if hosts:
        expected_results = [r for r in expected_results if r.get("host") in hosts]

    expected_map = {(r["case"], r["host"], r["frame"]): r["sha256"] for r in expected_results}
    actual_map = {(r["case"], r["host"], r["frame"]): r["sha256"] for r in actual["results"]}

    mismatches = []
    for key, expected_sha in expected_map.items():
        if key not in actual_map:
            mismatches.append(f"Missing case {key[0]} host {key[1]} frame {key[2]} in actual results")
        elif actual_map[key] != expected_sha:
            mismatches.append(
                f"Mismatch for case {key[0]} host {key[1]} frame {key[2]}: "
                f"expected {expected_sha}, got {actual_map[key]}"
            )
            
    if mismatches:
        for err in mismatches:
            print(err)
        raise AssertionError(f"Baseline mismatch in {path}")
    print("All verification cases successfully matched!")

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate or verify Neo-F3KDB frame hash baselines.")
    parser.add_argument("command", choices=["generate", "verify"])
    parser.add_argument("--plugin", required=True, type=Path)
    parser.add_argument("--avs-dump", type=Path)
    parser.add_argument("--cases-dir", type=Path, default=Path(__file__).parent / "cases")
    parser.add_argument("--golden", required=True, type=Path)
    parser.add_argument("--tier", default="smoke")
    parser.add_argument("--hosts", nargs="+", default=["vs", "avs"], choices=["vs", "avs"])
    parser.add_argument("--backend", default="purec", choices=["purec", "highway"])
    args = parser.parse_args(argv)

    plugin_path = args.plugin.resolve()
    avs_dump = args.avs_dump.resolve() if args.avs_dump else None
    cases = select_cases(load_cases(args.cases_dir), args.tier)
    results = run_cases(cases, args.hosts, plugin_path, avs_dump, args.backend)
    
    if args.command == "generate":
        write_golden(args.golden, results)
        print(f"Successfully generated golden baseline file: {args.golden}")
    else:
        verify_golden(args.golden, results, hosts=args.hosts)
    return 0

def run_cases(
    cases: list[dict],
    hosts: list[str],
    plugin_path: Path,
    avs_dump: Path | None,
    backend: str,
) -> list[dict]:
    results: list[dict] = []
    for case in cases:
        if "vs" in hosts:
            results.extend(run_vs_case(case, plugin_path, backend))
        if "avs" in hosts:
            if avs_dump is None:
                raise ValueError("--avs-dump is required when running AviSynth cases")
            results.extend(run_avs_case(case, plugin_path, avs_dump, backend))
    return results

if __name__ == "__main__":
    raise SystemExit(main())
