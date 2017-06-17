# $Id$
# Converts a log-file such that it can be used as reference in regression tests.
# Basically, it deletes lines which are pointless in comparison with references.

1,5 d
6 s/^ *[^ ][^ ]* //
7,/^ *VTF support/d
/^ *$/d
/^ *iter=/d
/^===============================================================/,$ d
/^====/d
/VTF/d
/umber of visualization points/d
/re-computing sparsity pattern/d
/^Writing/d
/^\.*$/d
/^ *Max/s/ node [0-9]*$//
# Since we use grep, some special characters that might occur must be escaped.
s/\*/\\\*/g
s/\[/\\\[/g
s/^ *--*//
/^Input file:/i\
