# gcFT8 (garbage collector for Franke-Taylor design, 8-FSK modulation in ham radio communication)
Automatic FT8 QSO systeme.
QRPLAB digital interface commande ligne FT8 tool.

QSO.log is the loogbook.
Call_Table.log for hashtable filter on callsign.

Mabe adaptation for ft817 in the futur.

Installation:

	git clone https://github.com/F4HTB/gcFT8
	make clean
	make

Example of use:
	clFT8 -d plughw:CARD=PCH,DEV=0 -C F4JJJ -L JN38 -F 14074000 -S /dev/ttyACM0 -x 1 -b

Command line option:

	-d sound device, prefer to use plughw
	-C your Callsign
	-L your locator
	-F TRX frequency
	-S serial device
	-x for set filter\n0 random (default)\n1 best decode score\n2 max distance\n3 min distance
	-b for console beep on log

	-x for set filter
		0 random (default)
		1 best decode score
		2 max distance
		3 min distance
		
	-b for console beep on log QSO after finish

	Red color: info
	Blue: CQ finded
	Magenta: filter by info missing or non standard message or call already made or Empty callsign
