# $Id$
# Converts a log-file such that it can be used as reference in regression tests.
# Basically, it deletes lines which are pointless in comparison with references.
# Remember to manually insert the command-line arguments as the first line.

1,/^ Executing command:/d
1,15 s/^ *[^ ][^ ]* //
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
