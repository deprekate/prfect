default:
	python3 -m pip install ../prfect/ --user

clean:
	rm -fr build
	rm -fr dist
	rm -fr prfect.egg-info
	python3 -m pip uninstall prfect -y

