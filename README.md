# gcFT8

Automatic FT8/FT4 command-line QSO tool for ham radio use.

## Logs

`gcFT8` writes QSOs in ADIF format to one file per local callsign:

```text
QSO_<CALLSIGN>.adif
```

Example:

```text
QSO_F4HTB.adif
```

At startup, the same ADIF file is loaded to build the in-memory callsign filter table:

- with `--band`, QSOs are filtered by `MODE` + `BAND`;
- with `--frequency`, QSOs are filtered by `MODE` + integer MHz from `FREQ`;
- the runtime hash table stores only remote callsigns.

## Build

```bash
make clean
make
```

## Example

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --beep
```

## Options

```text
--help
--mode <ft8|ft4>
--sound-device <device>
--callsign <callsign>
--locator <locator>
--frequency <hz>
--band <80|60|40|30|20|17|15|12|11|10>
--serial-device <device>
--filter <0..6>
--beep
```

`--frequency` and `--band` are mutually exclusive.

## Filters

```text
0  Listen only, no automatic TX
1  Random CQ selection
2  Best decode score
3  Maximum distance
4  Minimum distance
5  Maximum SNR
6  Minimum SNR
```

## Colors

```text
Red      Local station related message
Blue     CQ candidate
Magenta  Filtered CQ, missing info, non-standard message, already worked callsign or empty callsign
```

Before using it, check whether automated operation is allowed in your country. For testing only.
