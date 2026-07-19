if LoadPackage("qpa") <> true then
    Error("failed to load QPA");
fi;

q0 := Quiver(1, []);;
a0 := PathAlgebra(GF(2), q0);;
module := RightModuleOverPathAlgebra(a0, [1], []);;
lift := CALL_WITH_CATCH(LiftIdempotents, [module]);;

equality_q := Quiver(2, [[1, 2, "a"]]);;
equality_algebra := PathAlgebra(GF(2), equality_q);;
equal_left := RightModuleOverPathAlgebra(
    equality_algebra, [1, 1], [["a", [[Zero(GF(2))]]]]);;
equal_right := RightModuleOverPathAlgebra(
    equality_algebra, [1, 1], [["a", [[Zero(GF(2))]]]]);;
different_dimension := RightModuleOverPathAlgebra(
    equality_algebra, [2, 1],
    [["a", [[Zero(GF(2))], [Zero(GF(2))]]]]);;
different_matrix := RightModuleOverPathAlgebra(
    equality_algebra, [1, 1], [["a", [[One(GF(2))]]]]);;

q := Quiver(1, [[1, 1, "x"]]);;
parent := PathAlgebra(GF(2), q);;
ideal := Ideal(parent, [parent.x * parent.x]);;
gb := GroebnerBasis(ideal, [parent.x * parent.x]);;
algebra := parent / ideal;;
presentation := [[algebra.x]];;
resolution := ProjectiveResolutionFpPathAlgebraModule(
    algebra, ideal, presentation, 1);;
three_argument_resolution := CALL_WITH_CATCH(
    ProjectiveResolutionFpPathAlgebraModule,
    [algebra, ideal, presentation]);;
next_projective := CALL_WITH_CATCH(
    FindNextRProjective, [resolution, parent, 1]);;
next_syzygy := CALL_WITH_CATCH(FindNextSyzygy, [resolution, 1]);;

initial_r_vertices_bound := HasRProjectivesVertexList(resolution);;
initial_f_prime_bound := HasProjectivesFPrimeList(resolution);;
SetRProjectivesVertexList(resolution, [[q.v1]]);;
SetProjectivesFPrimeList(resolution, [[]]);;

doc_dirs := DirectoriesPackageLibrary("QPA", "doc");;
Print("{\"schema_version\":1,\"kind\":\"final_declared_surface\",\"oracle\":{");
Print("\"gap_version\":\"", GAPInfo.Version, "\"");
Print(",\"gap_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GAP_REVISION, "\"");
Print(",\"qpa_version\":\"", PackageInfo("qpa")[1].Version, "\"");
Print(",\"qpa_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_QPA_REVISION, "\"");
Print(",\"gbnp_version\":\"", PackageInfo("gbnp")[1].Version, "\"");
Print(",\"gbnp_revision\":\"",
    GAPInfo.SystemEnvironment.QPA_ORACLE_GBNP_REVISION, "\"}");
Print(",\"uninstalled_success\":[", lift[1], ",",
    three_argument_resolution[1], ",", next_projective[1], ",",
    next_syzygy[1], "]");
Print(",\"module_equality\":[", equal_left = equal_left, ",",
    IsIdenticalObj(equal_left, equal_left), ",", equal_left = equal_right,
    ",", IsIdenticalObj(equal_left, equal_right), ",",
    equal_left = different_dimension, ",", equal_left = different_matrix,
    ",", equal_left = module, "]");
Print(",\"initial_attribute_bound\":[", initial_r_vertices_bound, ",",
    initial_f_prime_bound, "]");
Print(",\"set_attribute_bound\":[",
    HasRProjectivesVertexList(resolution), ",",
    HasProjectivesFPrimeList(resolution), "]");
Print(",\"set_attribute_values\":[",
    RProjectivesVertexList(resolution) = [[q.v1]], ",",
    ProjectivesFPrimeList(resolution) = [[]], "]");
Print(",\"info\":[", IsBound(InfoProjectiveResolutionFpPathAlgebraModule),
    ",", InfoLevel(InfoProjectiveResolutionFpPathAlgebraModule), "]");
Print(",\"documentation\":{");
Print("\"function_bound\":", IsBound(MakeQPADocumentation));
Print(",\"directory_count\":", Length(doc_dirs));
Print(",\"main_exists\":",
    IsExistingFile(Filename(doc_dirs[1], "qpadocumentation.xml")));
Print(",\"main\":\"qpadocumentation.xml\",\"files\":[]");
Print(",\"book_name\":\"QPA\",\"renderer\":\"MathJax\"");
Print(",\"copy_html_style_files\":true}}");
Print("\n");

QUIT_GAP(0);
