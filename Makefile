all:
	@cd Hipo; make
	@cd Banks; make
	@cd Physics; make
clean:
	@cd Hipo; make clean
	@cd Banks; make clean
	@cd Physics; make clean
	@echo "Cleaning lib directory"
	@rm -rf lib/*
