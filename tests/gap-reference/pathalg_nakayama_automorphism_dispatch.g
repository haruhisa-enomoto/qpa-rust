if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

Succeeded := function(operation, arguments)
    return CALL_WITH_CATCH(operation, arguments)[1];
end;

cycle_parent := PathAlgebra(GF(2),
    Quiver(["c1", "c2", "c3"],
        [["c1", "c2", "x"], ["c2", "c3", "y"],
         ["c3", "c1", "z"]]));;
cycle := cycle_parent /
    [cycle_parent.x * cycle_parent.y,
     cycle_parent.y * cycle_parent.z,
     cycle_parent.z * cycle_parent.x];;
auto := NakayamaAutomorphism(cycle);;
form := FrobeniusForm(cycle);;
functional := FrobeniusLinearFunctional(cycle);;

foreign_parent := PathAlgebra(GF(2), Quiver(["f"], []));;
foreign := foreign_parent / [];;
foreign_element := One(foreign);;

a3_parent := PathAlgebra(GF(2),
    Quiver(["q1", "q2", "q3"],
        [["q1", "q2", "u"], ["q2", "q3", "v"]]));;
nonselfinjective := a3_parent / [a3_parent.u * a3_parent.v];;

nonadmissible_parent := PathAlgebra(GF(2),
    Quiver(["v1", "v2"], [["v1", "v2", "a"]]));;
nonadmissible := nonadmissible_parent / [nonadmissible_parent.a];;

Print("{\"schema_version\":1,\"source\":{",
      "\"gap\":\"4.17dev@4efe877859ab0bc0ac041bf635a996fb897ee0c9\",",
      "\"qpa\":\"1.37@8d5f847258d4853700ce3a5d947db2a1441109f1\"},",
      "\"foreign_auto_success\":", Succeeded(auto, [foreign_element]), ",",
      "\"foreign_form_left_success\":",
          Succeeded(form, [foreign_element, One(cycle)]), ",",
      "\"foreign_form_right_success\":",
          Succeeded(form, [One(cycle), foreign_element]), ",",
      "\"foreign_functional_success\":",
          Succeeded(functional, [foreign_element]), ",",
      "\"nonself_auto_is_false\":",
          NakayamaAutomorphism(nonselfinjective) = false, ",",
      "\"nonself_form_is_false\":",
          FrobeniusForm(nonselfinjective) = false, ",",
      "\"nonself_functional_is_false\":",
          FrobeniusLinearFunctional(nonselfinjective) = false, ",",
      "\"nonself_order_success\":",
          Succeeded(OrderOfNakayamaAutomorphism, [nonselfinjective]), ",",
      "\"nonadmissible_auto_success\":",
          Succeeded(NakayamaAutomorphism, [nonadmissible]), ",",
      "\"nonadmissible_form_success\":",
          Succeeded(FrobeniusForm, [nonadmissible]), ",",
      "\"nonadmissible_functional_success\":",
          Succeeded(FrobeniusLinearFunctional, [nonadmissible]), ",",
      "\"nonadmissible_order_success\":",
          Succeeded(OrderOfNakayamaAutomorphism, [nonadmissible]), "}\n");

QUIT_GAP(0);
