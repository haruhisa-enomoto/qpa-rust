if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

ar := ARQuiverNumerical(5, 3, [[0],[1,0],[4,0],[2,1],[3,4]]);;
PrintMultiplicityVector([0,2,0,1]);
PrintMultiplicityVector([0,0]);
PrintMultiplicityVectors([[1,0,3],[0,4,0]]);
Print(ar, "\n");
View(ar); Print("\n");

QUIT_GAP(0);
