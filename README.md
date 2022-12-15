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

    * 5: Similar to sample mode 4 but performing additional checks for better details preserving. (> r8)\
        For more info - https://forum.doom9.org/showthread.php?p=1652256#post1652256 (avgDif, maxDif, midDif1, midDif2)\
        `blur_first` doesn't have effect for this sample mode.

    Reference points are randomly picked within the `range`.

- *input_depth* (removed)

- *input_mode* (removed)

- *output_mode* (removed)

- *mt*

    Process planes in parallel. Default: true.

    If you notice a dead lock under extreme condition, try disabling it.

- *scale* (> r8)

    Whether to use threshold parameters (Y, Cb, Cr...) within the internal bit depth range (0..65535).
    
    Default: false.
    
- *Y_1 / Cb_1 / Cr_1 (maxDif),  Y_2 / Cb_2 / Cr_2 (midDif1, midDif2)* (> r8)

    Additional thresholds for `sample_mode=5`.    

## Compilation

```cmd
cmake -B build\x86 -S . -DCMAKE_GENERATOR_PLATFORM=Win32 -D_DIR=x86
cmake -B build\x64 -S . -DCMAKE_GENERATOR_PLATFORM=x64 -D_DIR=x64
cmake --build build\x86 --config Release
cmake --build build\x64 --config Release
```

## Compilation (GCC)

```bash
cmake -B build/gcc -S . -G "MSYS Makefiles" -D_DIR=gcc
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
