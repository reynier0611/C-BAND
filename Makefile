all:
	@cd Hipo; make
	@cd Banks; make
clean:
	@cd Hipo; make clean
	@cd Banks; make clean
	@echo "Cleaning lib directory"
	@rm -rf lib/*
