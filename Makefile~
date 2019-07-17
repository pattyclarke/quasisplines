
# The command "make preparsed" takes the sage files, preparses them, and tucks them into the package directory                                                 

preparsed:
	sage -preparse *.sage
	for file in *.sage.py; do mv "$$file" "$${file%.sage.py}.py"; done
	for file in *.py; do mv "$$file" "sagednd/$$file"; done


