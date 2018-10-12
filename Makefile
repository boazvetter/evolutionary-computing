.PHONY set-ld-path:
set-ld-path:
	export LD_LIBRARY_PATH=$$PWD

player66.class: contest.jar player66.java
	javac -cp contest.jar player66.java

submission: MainClass.txt $(wildcard *.class) player66.class
	jar cmf MainClass.txt submission.jar $(wildcard *.class)

.PHONY katsuura:
katsuura: submission testrun.jar set-ld-path
	java -jar testrun.jar -submission=player66 -evaluation=KatsuuraEvaluation -seed=1

.PHONY sphere:
sphere: submission testrun.jar
	java -jar testrun.jar -submission=player66 -evaluation=SphereEvaluation -seed=1

.PHONY bent-cigar:
bent-cigar: submission testrun.jar set-ld-path
	java -jar testrun.jar -submission=player66 -evaluation=BentCigarFunction -seed=1

.PHONY schaffers:
schaffers: submission testrun.jar
	java -jar testrun.jar -submission=player66 -evaluation=SchaffersEvaluation -seed=1


.PHONY test-all:
test-all: katsuura sphere bent-cigar schaffers

