function printcrop(f,opt,file)
% JPP 28.5.2020

print(f,opt,file)
unix(['pdfcrop ' file ' ' file]);

end

