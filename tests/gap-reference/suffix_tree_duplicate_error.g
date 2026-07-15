if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

tree := CreateSuffixTree();;
InsertPatternIntoSuffixTree(tree, [1,2,1]);;
InsertPatternIntoSuffixTree(tree, [2,1]);;
InsertPatternIntoSuffixTree(tree, [1,2]);;
InsertPatternIntoSuffixTree(tree, [3]);;
InsertPatternIntoSuffixTree(tree, [1,2,1]);;
DeletePatternFromSuffixTree(tree, 2);;
DeletePatternFromSuffixTree(tree, 1);;
AllPatternsInSequence(tree, [0,1,2,1,2,1,3]);;
QUIT_GAP(0);
