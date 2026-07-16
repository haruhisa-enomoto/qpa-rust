if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

table := EmptySCTable(3, 0);;
SetEntrySCTable(table, 1, 1, [1, 1]);;
SetEntrySCTable(table, 1, 2, [1, 2]);;
SetEntrySCTable(table, 2, 3, [1, 2]);;
SetEntrySCTable(table, 3, 3, [1, 3]);;
algebra := AlgebraByStructureConstants(Rationals, table);;
IsCommutative(algebra);;
opposite := OppositeAlgebra(algebra);;

Print(opposite, "\n");
QUIT_GAP(0);
