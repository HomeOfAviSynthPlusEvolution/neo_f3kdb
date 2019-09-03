# Neo f3kdb (forked from flash3kyuu_deband)

Neo f3kdb Copyright(C) 2019 msg7086

F3kdb is a deband filter originally written for AviUtl and later ported to AviSynth by [SAPikachu](https://github.com/SAPikachu) many years ago.

Legacy format support is removed and the filter is renamed to avoid confusion.

## Usage

```python
# AviSynth+
neo_f3kdb(clip, Y=64, Cb=64, Cr=64, grainY=0, grainC=0, mt=true)
```

[Check original usage documents.](https://f3kdb.readthedocs.io/en/stable/usage.html)


## Compilation

```cmd
mkdir build\x86
pushd build\x86
cmake -DCMAKE_GENERATOR_PLATFORM=Win32 -D_DIR=x86 ..\..\
popd
mkdir build\x64
pushd build\x64
cmake -DCMAKE_GENERATOR_PLATFORM=x64 -D_DIR=x64 ..\..\
popd
cmake --build build\x86 --config Release
cmake --build build\x64 --config Release
```

## Compilation (GCC)

```bash
mkdir -p build/gcc
pushd build/gcc
cmake -G "MSYS Makefiles" -D_DIR=gcc ../../
popd
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
