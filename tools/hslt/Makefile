all: deriv.txt

deriv.txt:a010.txt a100.txt b010.txt b100.txt c010.txt c100.txt
	cat a010.txt a100.txt b010.txt b100.txt c010.txt c100.txt > deriv.txt

MAT=a.max b.max c.max
FRAMES=010.max 100.max
a010.txt:vhslt.max $(MAT) $(FRAMES)
	cp a.max mat.max
	cp 010.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@
a100.txt:vhslt.max $(MAT) $(FRAMES)
	cp a.max mat.max
	cp 100.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@
b010.txt:vhslt.max $(MAT) $(FRAMES)
	cp b.max mat.max
	cp 010.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@
b100.txt:vhslt.max $(MAT) $(FRAMES)
	cp b.max mat.max
	cp 100.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@
c010.txt:vhslt.max $(MAT) $(FRAMES)
	cp c.max mat.max
	cp 010.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@
c100.txt:vhslt.max $(MAT) $(FRAMES)
	cp c.max mat.max
	cp 100.max frame.max
	maxima -q -b vhslt.max | grep -v "(%i" | sed 's/:/=/g' | sed 's/abs//g' | sed 's/block(.*],/const float /g' > $@

