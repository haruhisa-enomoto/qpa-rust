if LoadPackage("qpa") <> true then Error("failed to load QPA"); fi;

table := EmptySCTable(1, 0);;
SetEntrySCTable(table, 1, 1, [1, 1]);;
algebra := AlgebraByStructureConstants(Rationals, table);;
Print(OppositeAlgebra(algebra), "\n");
