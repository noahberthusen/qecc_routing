LoadPackage("QDistRnd");

filedir := DirectoriesPackageLibrary("QDistRnd","matrices");;
# output := OutputTextFile(Filename(filedir, "distances.txt"), false);;
# SetPrintFormattingStatus(output, false);

qx_filename := "QX.mtx";;
qz_filename := "QZ.mtx";;

lisX:=ReadMTXE(Filename(filedir, qx_filename), 0);;
lisZ:=ReadMTXE(Filename(filedir, qz_filename), 0);;
GX:=lisX[3];;
GZ:=lisZ[3];;
dz := DistRandCSS(GX,GZ,1000,0,0);;

# PrintTo(output, i);
# PrintTo(output, ",");
# PrintTo(output, dz);
# PrintTo(output, "\n");
Print(dz);