if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

table := EmptySCTable(0, 0);
algebra := AlgebraByStructureConstants(Rationals, table);
FindMultiplicativeIdentity(algebra);
