if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

tree := CreateSuffixTree();;
InsertPatternIntoSuffixTree(tree, []);;
QUIT_GAP(0);
