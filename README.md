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

`make` builds the release binary by default. Use `make debug` for an AddressSanitizer debug binary named `gcFT8-debug`.

## Example

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --max-same-tx-repeats 6 --beep
```

## Automatic CQ Selection Examples

Randomly select CQ calls at or above -18 dB:

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --snr-min -18
```

Select only CQ callsigns starting with `JA`, `VK`, or `ZL`:

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --only-prefix JA,VK,ZL
```

Select CQ calls in two-letter Maidenhead locator zones. `BP:CO` accepts `BO`, `BP`, `CO`, and `CP`:

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --only-locator-zone BP:CO
```

Select only special-purpose CQ tags such as `CQ POTA ...`, `CQ SOTA ...`, or `CQ DX ...`:

```bash
./gcFT8 --mode ft8 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 1 --only-sp-tag POTA,SOTA,DX
```

Combine filters before the final automatic CQ choice:

```bash
./gcFT8 --mode ft2 --sound-device plughw:CARD=PCH,DEV=0 --callsign F4JJJ --locator JN38 --band 20 --serial-device /dev/ttyACM0 --filter 5 --snr-min -16 --only-prefix JA,VK,ZL --only-sp-tag POTA,SOTA --only-locator-zone BP:FL,IO:KM,ON:PL
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
--max-same-tx-repeats <1..100>
--snr-min <snr>
--only-prefix <prefix[,prefix...]>
--only-sp-tag <tag[,tag...]>
--only-locator-zone <LL:LL[,LL:LL...]>
--beep
```

The following options are required: `--sound-device`, `--callsign`, `--locator`, `--serial-device`, and exactly one of `--frequency` or `--band`.

`--frequency` and `--band` are mutually exclusive. There is no default band or frequency.

`--max-same-tx-repeats` is optional and defaults to `6`. It limits repeated transmissions of the same QSO sequence, including repeated `seq 0` answers to a CQ from the same station. CQ messages blocked by this limit are displayed in gray.

`--snr-min` is optional. When present, automatic CQ selection ignores candidates below the given SNR, for example `--snr-min -18`.

`--only-prefix` is optional. When present, automatic CQ selection only keeps CQ callsigns matching one of the comma-separated prefixes, case-insensitively. Simple portable suffixes such as `/P`, `/M`, `/MM`, `/AM`, and `/QRP` are ignored for matching.

`--only-sp-tag` is optional. When present, automatic CQ selection only keeps tagged CQ calls matching one of the comma-separated special-purpose tags, for example `POTA`, `SOTA`, `IOTA`, `DX`, or `TEST`. Matching is exact after uppercase normalization; an untagged `CQ CALL LOC` is rejected by this filter.

`--only-locator-zone` is optional. When present, automatic CQ selection only keeps locators whose first two Maidenhead letters fall inside one of the inclusive ranges. For example, `BP:CO` accepts `BO`, `BP`, `CO`, and `CP`.

These optional filters only affect automatic CQ candidate selection. Decoded messages are still displayed, direct messages to your station are not blocked, and the ADIF already-worked filter still applies.

At startup, `gcFT8` prints a summary of the active mode, frequency, log file, filter mode, `--max-same-tx-repeats`, `--snr-min`, `--only-prefix`, `--only-sp-tag`, and `--only-locator-zone` settings.

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
Magenta  Filtered CQ, missing locator/callsign, or CQ rejected by optional filters
Gray     CQ rejected by max repeated TX attempts
```

Before using it, check whether automated operation is allowed in your country. For testing only.
