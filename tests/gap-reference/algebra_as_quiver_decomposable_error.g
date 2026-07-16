if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

field := GF(2);;
table := EmptySCTable(2, Zero(field));;
SetEntrySCTable(table, 1, 1, [One(field), 1, Zero(field), 2]);;
SetEntrySCTable(table, 1, 2, [Zero(field), 1, Zero(field), 2]);;
SetEntrySCTable(table, 2, 1, [Zero(field), 1, Zero(field), 2]);;
SetEntrySCTable(table, 2, 2, [Zero(field), 1, One(field), 2]);;
algebra := AlgebraByStructureConstants(field, table);;
AlgebraAsQuiverAlgebra(algebra);

QUIT_GAP(0);
