#!/bin/sh
# icover.sh NPTS  or  icover.sh i,j
if test ! -z "$2"
	then
	case "$1" in
	3)
		shift
		;;
	*)
		echo icover.sh 3 NPTS   or  icover.sh 3 i,j >&2
		exit 1
		;;
	esac
fi
N=$1
if test -z "$N"
	then echo icover.sh NPTS >&2
	exit 1
fi
sed  -n "
/^3 $N /,/^\$/p
" Sloane/codes.icover | grep -v ',' | Sloane/creconstruct
