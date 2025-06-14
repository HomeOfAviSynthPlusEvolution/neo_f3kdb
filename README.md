# Neo f3kdb (forked from flash3kyuu_deband)

Neo f3kdb Copyright(C) 2019-2020 Xinyue Lu, and previous developers

F3kdb is a deband filter. It was originally written for AviUtl by [bunyuchan](https://twitter.com/bunyuchan) and later ported to AviSynth by [SAPikachu](https://github.com/SAPikachu) many years ago.

Legacy format support was removed and a few options that are no longer useful were also removed. Due to API change, the project has been renamed from f3kdb to Neo_f3kdb to avoid confusion. SSE4.1 is now required to run optimized routine. SSE4.1 is supported since Intel Penryn (2007) and AMD bulldozer (2011). AVX routine didn't show much performance benefit and is not included.

## Usage

```python
# AviSynth+
LoadPlugin("neo-f3kdb.dll")
neo_f3kdb(clip, y=64, cb=64, cr=64, grainy=0, grainc=0, ...)
# VapourSynth
core.neo_f3kdb.Deband(clip, y=64, cb=64, cr=64, grainy=0, grainc=0, ...)
```

[Check original usage documents.](https://f3kdb.readthedocs.io/en/stable/usage.html)

- *sample_mode*

    * 1: Column references.

            +
            o
            +

    * 2: Square references.

            + +
             o
            + +

    * 3: Row references. (> r2)

            + o +

    * 4: Average of sample mode 1 and 3. (> r2)

             +
            (o) => A
             +

            + (o) + => B

            (A + B) / 2

    * 5: (Integer-based) Similar to sample mode 4 but uses multiple thresholds for detail preservation. (>r8)<br>
        Optimized for speed version of https://forum.doom9.org/showthread.php?p=1652256#post1652256.<br>
        `blur_first` doesn't have effect for this sample mode.<br>
        `Y`/`Cb`/`Cr` - for this mode they are used for the `avgDif` check – the difference between the current pixel and the average of all four cross-shaped reference pixels.

    * 6: (Floating-point) Similar to sample mode 4 but uses multiple thresholds for detail preservation. (>r9)<br>
        Direct implementation of https://forum.doom9.org/showthread.php?p=1652256#post1652256.<br>
        `blur_first` doesn't have effect for this sample mode.<br>
        `Y`/`Cb`/`Cr` - for this mode they are used for the `avgDif` check – the difference between the current pixel and the average of all four cross-shaped reference pixels.

    * 7: (Floating-point) An extension of sample_mode=6 that adds a gradient angle check for more intelligent detail preservation. (>r9)<br>
        Direct implementation of https://forum.doom9.org/showthread.php?p=1652256#post1652256.<br>
        `blur_first` doesn't have effect for this sample mode.<br>
        `Y`/`Cb`/`Cr` - for this mode they are used for the `avgDif` check – the difference between the current pixel and the average of all four cross-shaped reference pixels.

    Reference points are randomly picked within the `range`.

- *input_depth* (removed)

- *input_mode* (removed)

- *output_mode* (removed)

- *opt*

    Sets which cpu optimizations to use.

    `sample_mode=1`, `sample_mode=2`, `sample_mode=3`, and `sample_mode=4` have `C++` and `SSE4.1` code.

    `sample_mode=5`, `sample_mode=6` and `sample_mode=7` have `C++`, `SSE4.1`, `AVX2` and `AVX-512` code.

    - `-1`: Auto-detect.
    - `0`: Use C++ code.
    - `1`: Use SSE4.1 code.
    - `2`: Use AVX2 code.
    - `3`: Use AVX-512 code.

    Default: `-1`.

- *mt*

    Process planes in parallel. Default: true.

    If you notice a dead lock under extreme condition, try disabling it.

- *scale* (> r8)

    Whether to use threshold parameters (Y, Cb, Cr...) within the internal bit depth range (0..65535).

    Default: false.

- *Y_1 / Cb_1 / Cr_1 (maxDif)* (> r8)

    Detail protection threshold (max difference) for `sample_mode=5`, `sample_mode=6` and `sample_mode=7`.

    This threshold applies to the `maxDif` check. `maxDif` is the largest absolute difference found between the current pixel and any of its four individual cross-shaped reference pixels. If this `maxDif` is greater than or equal to `Y_1`/`Cb_1`/`Cr_1`, the pixel is considered detail.

    Helps protect sharp edges and fine details from being blurred by the debanding process.

    The valid range is same as `Y`/`Cb`/`Cr`.

    Default value - they are equal to `Y`/`Cb`/`Cr`.

- *Y_2 / Cb_2 / Cr_2 (midDifs)* (> r8)

    Gradient/Texture protection threshold (mid-pair difference) for `sample_mode=5`, `sample_mode=6` and `sample_mode=7`.

    This threshold applies to the `midDif` checks. `midDif` measures how much the current pixel deviates from the midpoint of a pair of opposing reference pixels (one check for the vertical pair, one for the horizontal pair). If the current pixel is far from this midpoint (i.e., `midDif` is greater than or equal to `Y_2` / `Cb_2` / `Cr_2`), it might indicate a texture.

    This helps distinguish true banding in gradients from textured areas or complex details.

    The valid range is same as `Y`/`Cb`/`Cr`.

    Default value - they are equal to `Y`/`Cb`/`Cr`.

- *angle_boost* (>r9)

    A multiplier used in `sample_mode=7` to increase the debanding strength on consistent gradients.

    When the gradient angle check passes, the `Y`/`Cb`/`Cr`, `Y_1`/`Cb_1`/`Cr_1`, and `Y_2`/`Cb_2`/`Cr_2` thresholds are multiplied by this factor.

    A value greater than `1.0` boosts the strength. A value of `1.0` has no effect.

    Must be a positive number.

    Default value - `1.5`.

- *max_angle* (>r9)

    The threshold for the gradient angle check in `sample_mode=7`.

    It represents the maximum allowed difference between the gradient angle of the center pixel and its reference pixels for the `angle_boost` to be applied. The gradient angle is normalized to a `[0.0, 1.0]` range.

    A smaller value is stricter and requires a more consistent gradient. A larger value is more lenient.

    The valid range is `0.0` to `1.0`.

    Default value - `0.15`.

## Compilation

```cmd
cmake -B build\x86 -S . -DCMAKE_GENERATOR_PLATFORM=Win32 -D_DIR=x86
cmake -B build\x64 -S . -DCMAKE_GENERATOR_PLATFORM=x64 -D_DIR=x64
cmake --build build\x86 --config Release
cmake --build build\x64 --config Release
```

## Compilation (GCC, Windows)

```bash
cmake -B build/gcc -S . -G "MSYS Makefiles" -D_DIR=gcc
cmake --build build/gcc
```

## Compilation (GCC, Unix-like)

```bash
cmake -B build/gcc -S . -G "Unix Makefiles" -D_DIR=gcc
cmake --build build/gcc
```

## License

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
