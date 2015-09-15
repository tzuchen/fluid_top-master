# type "make" command in Unix to create asme2e.pdf file 
all:
	latex dsccpaper
	bibtex dsccpaper
	latex dsccpaper
	latex dsccpaper
	dvips -Ppdf -G0 -tletter dsccpaper
	ps2pdf -dCompatibilityLevel=1.4 -dMAxSubsetPct=100 -dSubsetFonts=true -dEmbedAllFonts=true -sPAPERSIZE=letter dsccpaper.ps
#	dvips -o dsccpaper.ps dsccpaper
#	ps2pdf dsccpaper.ps dsccpaper.pdf
#	evince dsccpaper.pdf

abstract:
	latex dsccabstract
	dvips -o dsccabstract.ps dsccabstract
	ps2pdf dsccabstract.ps dsccabstract.pdf

clean:
	(rm -rf *.ps *.log *.dvi *.aux *.*% *.lof *.lop *.lot *.toc *.idx *.ilg *.ind *.bbl *blg)
