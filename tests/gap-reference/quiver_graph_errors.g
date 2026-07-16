if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

quiver := Quiver(1, []);;
foreign := Quiver(["x"], []);;
nonvertex := CALL_WITH_CATCH(FullSubquiver, [quiver, [1]]);;
foreign_vertex := CALL_WITH_CATCH(FullSubquiver, [quiver, [foreign.x]]);;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"nonvertex_success\":", nonvertex[1], ",",
      "\"foreign_vertex_success\":", foreign_vertex[1], "}\n");

QUIT_GAP(0);
