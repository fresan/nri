
all: example

example: randomindex.cc main.cc randomindex.h
	g++ -O2 -o example randomindex.cc main.cc
	
clean:
	rm -rf example

# EOF
