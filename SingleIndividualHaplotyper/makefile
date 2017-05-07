all: clean compile jar

clean: 
	rm -f SIH.jar
	rm -rf bin
	
compile:
	mkdir bin 
	javac -d bin src/mpg/molgen/sih/*/*.java

jar: 
	jar -cvfe SIH.jar mpg.molgen.sih.main.SIH -C bin mpg

