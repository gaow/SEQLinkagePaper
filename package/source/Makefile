pyinstaller=/home/usr/pyinstaller/pyinstaller.py
dist=/home/usr/bin

help:
	@echo " "
	@echo "Type ...............	To ................"
	@echo "make source		get source code release"
	@echo "make bin		get a standalone binary release"
	@echo "make clean		to clean up"
	@echo " "

source:
	cat META libslinco.py slinco.py | awk 'NR <=3 {print $$0; next}; {if ($$1 != "#" && $$NF != "#") print $$0}' > slinco
	chmod +x slinco
	mv slinco $(dist)

bin:
	cat META libslinco.py slinco.py | awk 'NR <=3 {print $$0; next}; {if ($$1 != "#" && $$NF != "#") print $$0}' > slinco
	$(pyinstaller) -F --distpath . slinco
	zip slinco-`grep "^VERSION" libslinco.py | sed -r 's/[^\"]*[\"]([^\"]*)[\"]/\1/g'`-Linux-`uname -m`.zip slinco
	make clean

clean:
	rm -rf slinco slinco.spec build dist
