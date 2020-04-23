
TESTS=$(shell find . -maxdepth 1 -type d -name '????_*' -printf '%p\n')
TEST_RESULTS=$(shell find . -maxdepth 1 -type d -name '????_*' -printf '%p.res\n')

.PHONY:$(TESTS) all clean clean_timings

all: report.html

clean:
	rm -f *.res

report.html: header.html footer.html $(TEST_RESULTS)
	cat header.html *.res footer.html > report.html

# create results for individual test:
%.res:
	@echo "making $*"
	@./createres.sh $*

# this rule forces rebuilding of one individual regression test.
# you should pass it the number only, i.e. `make 0051'.
%:
	@echo "forcing make of test $*"
	rm -f $*_*.res
	make $*_*.res
	make report.html

clean_timings:
	for i in $(TESTS); do rm -f $$i/wallclock; done
	for i in $(TESTS); do rm -f $$i/usertime; done
