# gcFT8

Automatic FT8/FT4/FT2 command-line QSO tool for ham radio use.

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

## Automatic CQ Selection Examples

Randomly select CQ calls at or above -18 dB:

```bash
./gcFT8 --mode ft8 --band 20 --filter 1 --snr-min -18
```

Select only CQ callsigns starting with `JA`, `VK`, or `ZL`:

```bash
./gcFT8 --mode ft8 --band 20 --filter 1 --only-prefix JA,VK,ZL
```

Select CQ calls in two-letter Maidenhead locator zones. `BP:CO` accepts `BO`, `BP`, `CO`, and `CP`:

```bash
./gcFT8 --mode ft8 --band 20 --filter 1 --only-locator-zone BP:CO
```

Combine filters before the final automatic CQ choice:

```bash
./gcFT8 --mode ft2 --band 20 --filter 5 --snr-min -16 --only-prefix JA,VK,ZL --only-locator-zone BP:FL,IO:KM,ON:PL
```

## Options

```text
--help
--mode <ft8|ft4|ft2>
--sound-device <device>
--callsign <callsign>
--locator <locator>
--frequency <hz>
--band <80|60|40|30|20|17|15|12|11|10>
--serial-device <device>
--filter <0..6>
--snr-min <snr>
--only-prefix <prefix[,prefix...]>
--only-locator-zone <LL:LL[,LL:LL...]>
--beep
```

`--frequency` and `--band` are mutually exclusive.

`--snr-min` is optional. When present, automatic CQ selection ignores candidates below the given SNR, for example `--snr-min -18`.

`--only-prefix` is optional. When present, automatic CQ selection only keeps CQ callsigns matching one of the comma-separated prefixes, case-insensitively. Simple portable suffixes such as `/P`, `/M`, `/MM`, `/AM`, and `/QRP` are ignored for matching.

`--only-locator-zone` is optional. When present, automatic CQ selection only keeps locators whose first two Maidenhead letters fall inside one of the inclusive ranges. For example, `BP:CO` accepts `BO`, `BP`, `CO`, and `CP`.

These optional filters only affect automatic CQ candidate selection. Decoded messages are still displayed, direct messages to your station are not blocked, and the ADIF already-worked filter still applies.

At startup, `gcFT8` prints a summary of the active mode, frequency, log file, filter mode, `--snr-min`, `--only-prefix`, and `--only-locator-zone` settings.

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
Cyan     RX/TX slot separator
Red      Local station related message
Blue     CQ candidate
Yellow   Already worked callsign
Magenta  Filtered CQ, missing info, non-standard message or empty callsign
```

Before using it, check whether automated operation is allowed in your country. For testing only.
