Search.setIndex({docnames:["classes/spatialpy","classes/spatialpy.core","classes/spatialpy.solvers","classes/spatialpy.stochss","index"],envversion:{"sphinx.domains.c":2,"sphinx.domains.changeset":1,"sphinx.domains.citation":1,"sphinx.domains.cpp":4,"sphinx.domains.index":1,"sphinx.domains.javascript":2,"sphinx.domains.math":2,"sphinx.domains.python":3,"sphinx.domains.rst":2,"sphinx.domains.std":2,"sphinx.ext.todo":2,"sphinx.ext.viewcode":1,sphinx:56},filenames:["classes/spatialpy.rst","classes/spatialpy.core.rst","classes/spatialpy.solvers.rst","classes/spatialpy.stochss.rst","index.rst"],objects:{"":[[0,0,0,"-","spatialpy"]],"spatialpy.core":[[1,1,1,"","BoundaryCondition"],[1,3,1,"","BoundaryConditionError"],[1,1,1,"","DataFunction"],[1,3,1,"","DataFunctionError"],[1,1,1,"","Domain"],[1,3,1,"","DomainError"],[1,1,1,"","Geometry"],[1,1,1,"","GeometryAll"],[1,3,1,"","GeometryError"],[1,1,1,"","GeometryExterior"],[1,1,1,"","GeometryInterior"],[1,1,1,"","InitialCondition"],[1,3,1,"","InitialConditionError"],[1,1,1,"","KDTree"],[1,1,1,"","Model"],[1,3,1,"","ModelError"],[1,1,1,"","OrderedDict"],[1,1,1,"","Parameter"],[1,3,1,"","ParameterError"],[1,1,1,"","PlaceInitialCondition"],[1,1,1,"","Reaction"],[1,3,1,"","ReactionError"],[1,1,1,"","Result"],[1,3,1,"","ResultError"],[1,1,1,"","ScatterInitialCondition"],[1,3,1,"","SimulationError"],[1,3,1,"","SimulationTimeout"],[1,1,1,"","Species"],[1,3,1,"","SpeciesError"],[1,1,1,"","UniformInitialCondition"],[1,1,1,"","VTKReader"],[1,3,1,"","VTKReaderError"],[1,3,1,"","VTKReaderIOError"],[1,0,0,"-","boundarycondition"],[1,0,0,"-","cleanup"],[1,6,1,"","cleanup_build_files"],[1,6,1,"","cleanup_core_files"],[1,6,1,"","cleanup_result_files"],[1,6,1,"","cleanup_tempfiles"],[1,1,1,"","combinations"],[1,0,0,"-","datafunction"],[1,0,0,"-","domain"],[1,0,0,"-","geometry"],[1,6,1,"","init_notebook_mode"],[1,0,0,"-","initialcondition"],[1,6,1,"","iplot"],[1,0,0,"-","model"],[1,0,0,"-","parameter"],[1,0,0,"-","reaction"],[1,0,0,"-","result"],[1,0,0,"-","spatialpyerror"],[1,0,0,"-","species"],[1,0,0,"-","vtkreader"]],"spatialpy.core.BoundaryCondition":[[1,2,1,"","expression"]],"spatialpy.core.DataFunction":[[1,2,1,"","map"]],"spatialpy.core.Domain":[[1,2,1,"","add_point"],[1,2,1,"","calculate_vol"],[1,2,1,"","closest_vertex"],[1,2,1,"","compile_prep"],[1,2,1,"","coordinates"],[1,2,1,"","create_2D_domain"],[1,2,1,"","create_3D_domain"],[1,2,1,"","distance_between_2_vertices"],[1,2,1,"","fill_with_particles"],[1,2,1,"","find_boundary_points"],[1,2,1,"","find_h"],[1,2,1,"","get_bounding_box"],[1,2,1,"","get_domain_size"],[1,2,1,"","get_num_voxels"],[1,2,1,"","get_type_def"],[1,2,1,"","get_vol"],[1,2,1,"","import_meshio_object"],[1,2,1,"","plot_types"],[1,2,1,"","read_msh_file"],[1,2,1,"","read_stochss_domain"],[1,2,1,"","read_stochss_subdomain_file"],[1,2,1,"","read_xml_mesh"],[1,2,1,"","set_properties"]],"spatialpy.core.Geometry":[[1,2,1,"","inside"]],"spatialpy.core.GeometryAll":[[1,2,1,"","inside"]],"spatialpy.core.GeometryExterior":[[1,2,1,"","inside"]],"spatialpy.core.GeometryInterior":[[1,2,1,"","inside"]],"spatialpy.core.InitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.KDTree":[[1,2,1,"","count_neighbors"],[1,4,1,"","data"],[1,1,1,"","innernode"],[1,1,1,"","leafnode"],[1,4,1,"","leafsize"],[1,4,1,"","m"],[1,4,1,"","maxes"],[1,4,1,"","mins"],[1,4,1,"","n"],[1,1,1,"","node"],[1,2,1,"","query"],[1,2,1,"","query_ball_point"],[1,2,1,"","query_ball_tree"],[1,2,1,"","query_pairs"],[1,4,1,"","size"],[1,2,1,"","sparse_distance_matrix"],[1,5,1,"","tree"]],"spatialpy.core.KDTree.innernode":[[1,5,1,"","children"],[1,5,1,"","split"],[1,5,1,"","split_dim"]],"spatialpy.core.KDTree.leafnode":[[1,5,1,"","children"],[1,5,1,"","idx"]],"spatialpy.core.Model":[[1,2,1,"","add_boundary_condition"],[1,2,1,"","add_data_function"],[1,2,1,"","add_domain"],[1,2,1,"","add_initial_condition"],[1,2,1,"","add_parameter"],[1,2,1,"","add_reaction"],[1,2,1,"","add_species"],[1,2,1,"","compile_prep"],[1,2,1,"","delete_all_parameters"],[1,2,1,"","delete_all_reactions"],[1,2,1,"","delete_all_species"],[1,2,1,"","delete_parameter"],[1,2,1,"","delete_reaction"],[1,2,1,"","delete_species"],[1,2,1,"","get_all_parameters"],[1,2,1,"","get_all_reactions"],[1,2,1,"","get_all_species"],[1,2,1,"","get_num_reactions"],[1,2,1,"","get_num_species"],[1,2,1,"","get_parameter"],[1,2,1,"","get_reaction"],[1,2,1,"","get_species"],[1,4,1,"","reserved_names"],[1,2,1,"","run"],[1,2,1,"","sanitized_data_function_names"],[1,2,1,"","sanitized_parameter_names"],[1,2,1,"","sanitized_species_names"],[1,2,1,"","set_timesteps"],[1,4,1,"","special_characters"],[1,2,1,"","timespan"],[1,2,1,"","update_namespace"]],"spatialpy.core.OrderedDict":[[1,2,1,"","clear"],[1,2,1,"","copy"],[1,2,1,"","fromkeys"],[1,2,1,"","items"],[1,2,1,"","keys"],[1,2,1,"","move_to_end"],[1,2,1,"","pop"],[1,2,1,"","popitem"],[1,2,1,"","setdefault"],[1,2,1,"","update"],[1,2,1,"","values"]],"spatialpy.core.PlaceInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.Reaction":[[1,2,1,"","add_product"],[1,2,1,"","add_reactant"],[1,2,1,"","annotate"],[1,2,1,"","initialize"]],"spatialpy.core.Result":[[1,2,1,"","export_to_csv"],[1,2,1,"","get_property"],[1,2,1,"","get_species"],[1,2,1,"","get_timespan"],[1,2,1,"","plot_property"],[1,2,1,"","plot_species"],[1,2,1,"","read_step"]],"spatialpy.core.ScatterInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.Species":[[1,2,1,"","set_diffusion_coefficient"]],"spatialpy.core.UniformInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.VTKReader":[[1,2,1,"","get_array_name"],[1,2,1,"","get_arrays"],[1,2,1,"","get_num_points"],[1,2,1,"","get_points"],[1,2,1,"","read_file"],[1,2,1,"","set_filename"]],"spatialpy.core.boundarycondition":[[1,1,1,"","BoundaryCondition"]],"spatialpy.core.boundarycondition.BoundaryCondition":[[1,2,1,"","expression"]],"spatialpy.core.cleanup":[[1,6,1,"","cleanup_build_files"],[1,6,1,"","cleanup_core_files"],[1,6,1,"","cleanup_result_files"],[1,6,1,"","cleanup_tempfiles"]],"spatialpy.core.datafunction":[[1,1,1,"","DataFunction"]],"spatialpy.core.datafunction.DataFunction":[[1,2,1,"","map"]],"spatialpy.core.domain":[[1,1,1,"","Domain"]],"spatialpy.core.domain.Domain":[[1,2,1,"","add_point"],[1,2,1,"","calculate_vol"],[1,2,1,"","closest_vertex"],[1,2,1,"","compile_prep"],[1,2,1,"","coordinates"],[1,2,1,"","create_2D_domain"],[1,2,1,"","create_3D_domain"],[1,2,1,"","distance_between_2_vertices"],[1,2,1,"","fill_with_particles"],[1,2,1,"","find_boundary_points"],[1,2,1,"","find_h"],[1,2,1,"","get_bounding_box"],[1,2,1,"","get_domain_size"],[1,2,1,"","get_num_voxels"],[1,2,1,"","get_type_def"],[1,2,1,"","get_vol"],[1,2,1,"","import_meshio_object"],[1,2,1,"","plot_types"],[1,2,1,"","read_msh_file"],[1,2,1,"","read_stochss_domain"],[1,2,1,"","read_stochss_subdomain_file"],[1,2,1,"","read_xml_mesh"],[1,2,1,"","set_properties"]],"spatialpy.core.geometry":[[1,1,1,"","Geometry"],[1,1,1,"","GeometryAll"],[1,1,1,"","GeometryExterior"],[1,1,1,"","GeometryInterior"]],"spatialpy.core.geometry.Geometry":[[1,2,1,"","inside"]],"spatialpy.core.geometry.GeometryAll":[[1,2,1,"","inside"]],"spatialpy.core.geometry.GeometryExterior":[[1,2,1,"","inside"]],"spatialpy.core.geometry.GeometryInterior":[[1,2,1,"","inside"]],"spatialpy.core.initialcondition":[[1,1,1,"","InitialCondition"],[1,1,1,"","PlaceInitialCondition"],[1,1,1,"","ScatterInitialCondition"],[1,1,1,"","UniformInitialCondition"]],"spatialpy.core.initialcondition.InitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.initialcondition.PlaceInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.initialcondition.ScatterInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.initialcondition.UniformInitialCondition":[[1,2,1,"","apply"]],"spatialpy.core.model":[[1,1,1,"","Model"],[1,6,1,"","export_StochSS"]],"spatialpy.core.model.Model":[[1,2,1,"","add_boundary_condition"],[1,2,1,"","add_data_function"],[1,2,1,"","add_domain"],[1,2,1,"","add_initial_condition"],[1,2,1,"","add_parameter"],[1,2,1,"","add_reaction"],[1,2,1,"","add_species"],[1,2,1,"","compile_prep"],[1,2,1,"","delete_all_parameters"],[1,2,1,"","delete_all_reactions"],[1,2,1,"","delete_all_species"],[1,2,1,"","delete_parameter"],[1,2,1,"","delete_reaction"],[1,2,1,"","delete_species"],[1,2,1,"","get_all_parameters"],[1,2,1,"","get_all_reactions"],[1,2,1,"","get_all_species"],[1,2,1,"","get_num_reactions"],[1,2,1,"","get_num_species"],[1,2,1,"","get_parameter"],[1,2,1,"","get_reaction"],[1,2,1,"","get_species"],[1,4,1,"","reserved_names"],[1,2,1,"","run"],[1,2,1,"","sanitized_data_function_names"],[1,2,1,"","sanitized_parameter_names"],[1,2,1,"","sanitized_species_names"],[1,2,1,"","set_timesteps"],[1,4,1,"","special_characters"],[1,2,1,"","timespan"],[1,2,1,"","update_namespace"]],"spatialpy.core.parameter":[[1,1,1,"","Parameter"]],"spatialpy.core.reaction":[[1,1,1,"","Reaction"]],"spatialpy.core.reaction.Reaction":[[1,2,1,"","add_product"],[1,2,1,"","add_reactant"],[1,2,1,"","annotate"],[1,2,1,"","initialize"]],"spatialpy.core.result":[[1,1,1,"","Result"]],"spatialpy.core.result.Result":[[1,2,1,"","export_to_csv"],[1,2,1,"","get_property"],[1,2,1,"","get_species"],[1,2,1,"","get_timespan"],[1,2,1,"","plot_property"],[1,2,1,"","plot_species"],[1,2,1,"","read_step"]],"spatialpy.core.spatialpyerror":[[1,3,1,"","BoundaryConditionError"],[1,3,1,"","DataFunctionError"],[1,3,1,"","DomainError"],[1,3,1,"","GeometryError"],[1,3,1,"","InitialConditionError"],[1,3,1,"","ModelError"],[1,3,1,"","ParameterError"],[1,3,1,"","ReactionError"],[1,3,1,"","ResultError"],[1,3,1,"","SimulationError"],[1,3,1,"","SimulationTimeout"],[1,3,1,"","SpeciesError"],[1,3,1,"","VTKReaderError"],[1,3,1,"","VTKReaderIOError"]],"spatialpy.core.species":[[1,1,1,"","Species"]],"spatialpy.core.species.Species":[[1,2,1,"","set_diffusion_coefficient"]],"spatialpy.core.vtkreader":[[1,1,1,"","VTKReader"]],"spatialpy.core.vtkreader.VTKReader":[[1,2,1,"","get_array_name"],[1,2,1,"","get_arrays"],[1,2,1,"","get_num_points"],[1,2,1,"","get_points"],[1,2,1,"","read_file"],[1,2,1,"","set_filename"]],"spatialpy.solvers":[[2,0,0,"-","build_expression"],[2,0,0,"-","solver"]],"spatialpy.solvers.build_expression":[[2,1,1,"","BuildExpression"],[2,1,1,"","CppConverter"],[2,1,1,"","ExpressionConverter"],[2,1,1,"","ExpressionResults"],[2,1,1,"","PythonConverter"]],"spatialpy.solvers.build_expression.BuildExpression":[[2,1,1,"","ValidationVisitor"],[2,2,1,"","getexpr_cpp"],[2,2,1,"","getexpr_python"],[2,2,1,"","map_operator"],[2,4,1,"","operator_map"],[2,2,1,"","validate"],[2,2,1,"","with_blacklist"],[2,2,1,"","with_namespace"]],"spatialpy.solvers.build_expression.BuildExpression.ValidationVisitor":[[2,2,1,"","visit_Assign"],[2,2,1,"","visit_BinOp"],[2,2,1,"","visit_BoolOp"],[2,2,1,"","visit_Call"],[2,2,1,"","visit_Compare"],[2,2,1,"","visit_Name"],[2,2,1,"","visit_UnaryOp"]],"spatialpy.solvers.build_expression.CppConverter":[[2,1,1,"","CppExpressionTransformer"],[2,2,1,"","get_str"],[2,2,1,"","visit_And"],[2,2,1,"","visit_Or"]],"spatialpy.solvers.build_expression.CppConverter.CppExpressionTransformer":[[2,2,1,"","visit_BinOp"]],"spatialpy.solvers.build_expression.ExpressionConverter":[[2,2,1,"","convert_str"],[2,2,1,"","get_str"],[2,2,1,"","parse_comparison"],[2,2,1,"","parse_logical"],[2,2,1,"","parse_operator"],[2,2,1,"","visit_Add"],[2,2,1,"","visit_BinOp"],[2,2,1,"","visit_BoolOp"],[2,2,1,"","visit_Bytes"],[2,2,1,"","visit_Call"],[2,2,1,"","visit_Compare"],[2,2,1,"","visit_Constant"],[2,2,1,"","visit_Div"],[2,2,1,"","visit_Ellipsis"],[2,2,1,"","visit_Eq"],[2,2,1,"","visit_Gt"],[2,2,1,"","visit_GtE"],[2,2,1,"","visit_Lt"],[2,2,1,"","visit_LtE"],[2,2,1,"","visit_Mult"],[2,2,1,"","visit_Name"],[2,2,1,"","visit_NameConstant"],[2,2,1,"","visit_NotEq"],[2,2,1,"","visit_Num"],[2,2,1,"","visit_Pow"],[2,2,1,"","visit_Str"],[2,2,1,"","visit_Sub"],[2,2,1,"","visit_USub"],[2,2,1,"","visit_UnaryOp"]],"spatialpy.solvers.build_expression.PythonConverter":[[2,2,1,"","visit_And"],[2,2,1,"","visit_Or"]],"spatialpy.solvers.solver":[[2,1,1,"","Solver"]],"spatialpy.solvers.solver.Solver":[[2,2,1,"","compile"],[2,2,1,"","run"]],"spatialpy.stochss":[[3,0,0,"-","stochss_export"]],"spatialpy.stochss.stochss_export":[[3,6,1,"","export"]],spatialpy:[[1,0,0,"-","core"],[2,0,0,"-","solvers"],[3,0,0,"-","stochss"]]},objnames:{"0":["py","module","Python module"],"1":["py","class","Python class"],"2":["py","method","Python method"],"3":["py","exception","Python exception"],"4":["py","attribute","Python attribute"],"5":["py","property","Python property"],"6":["py","function","Python function"]},objtypes:{"0":"py:module","1":"py:class","2":"py:method","3":"py:exception","4":"py:attribute","5":"py:property","6":"py:function"},terms:{"0":[1,2],"0012333":1,"04":1,"0511773":1,"1":[1,2],"10":1,"11":1,"12295571":1,"13":1,"14":1,"15":1,"168474926":1,"19":1,"1993":1,"1993apj":1,"1999":1,"1d":1,"2":1,"20":1,"2000":1,"2002":1,"2003mnra":1,"2005":1,"2015":1,"2019":[0,1,2,3],"2022":[0,1,2,3],"2236068":1,"2333084":1,"23606798":1,"24617575":1,"26836782":1,"28942611":1,"29571802":1,"2d":1,"2df":1,"2r01eb014877":4,"3":[0,1,2,3,4],"300":1,"30088013":1,"31994999":1,"32580855":1,"33630734":1,"346":1,"37327919":1,"39923701":1,"3d":[1,4],"4":1,"412":1,"5":1,"500":1,"55396483":1,"5647":1,"56906522":1,"57714465":1,"6":1,"600":1,"6395831":1,"6473269":1,"64l":1,"7":1,"71124834":1,"7225693":1,"72658602":1,"75679479":1,"78h":1,"79428925":1,"8":1,"800":1,"80622577":1,"8658745":1,"87665969":1,"9":1,"90785663":1,"abstract":[1,2],"byte":2,"case":1,"class":[1,2],"default":1,"do":1,"export":[1,3],"float":1,"function":[1,2],"import":1,"int":[1,2],"long":1,"new":[1,2],"public":[0,1,2,3,4],"return":[1,2,3],"short":4,"static":[1,4],"super":2,"switch":1,"true":[1,2],A:[0,1,2,3,4],And:2,As:1,By:[1,2],FOR:[0,1,2,3],For:[1,4],If:[0,1,2,3,4],In:1,It:1,No:1,Not:2,Or:2,The:[1,2,4],To:1,Will:1,ab:1,abl:1,about:1,abov:1,absolut:1,acceler:1,accept:2,across:1,action:1,actual:1,ad:[1,2],adapt:1,add:[1,2],add_boundary_condit:1,add_data_funct:1,add_domain:1,add_initial_condit:1,add_paramet:1,add_point:1,add_product:1,add_react:1,add_speci:1,addit:[1,2],addition:[1,2],addtion:[1,2],adsab:1,advect:[0,1,2,3,4],against:[1,2],al:1,algorithm:[1,4],align:1,all:[1,2],allow:2,along:[0,1,2,3],alreadi:1,also:[1,2,4],alwai:[1,2],amount:1,an:[1,2,4],analysi:1,angular:1,ani:[0,1,2,3,4],anim:1,animation_attribut:1,animation_opt:1,annot:1,append:2,appli:1,approxim:1,ar:[1,2],argument:1,around:1,arrai:1,arraid:1,array_lik:1,articl:4,artifici:1,arxiv:1,asarrai:1,ascii:1,assembl:1,assign:[1,2],associ:1,assum:2,ast:2,astro:1,astronom:1,astrophys:1,attempt:1,attribut:1,author:4,auto_plai:1,automat:[1,2],autorang:1,avail:[1,2,4],averag:1,axi:1,background:1,balanced_tre:1,bare:2,base:[1,2,4],basic:1,bass:1,becaus:1,becom:1,been:1,befor:1,begin:1,behavior:[1,2],below:1,between:[1,2],bia:1,bin:1,binari:1,binboolopop:2,binop:2,biochem:1,bitand:2,bitor:2,bitxor:2,blacklist:2,blob:1,blue:1,bodi:1,bogu:1,bool:[1,2,3],boolop:2,both:4,bottom:1,bound:1,bound_cond:1,boundari:1,boundarycondit:0,boundaryconditionerror:1,box:1,boxsiz:1,bracket:1,branch:1,brian:4,broadcast:1,browser:1,brute:1,bugn:1,build:[1,2,4],build_dir:1,build_express:0,buildexpress:2,built:[1,2],bulk:1,bupu:1,bvf_phi:1,c0:1,c:[0,1,2,3],c_:1,calcul:1,calculate_vol:1,call:[1,2],can:[0,1,2,3,4],caus:[1,2],cdn:1,cell:1,charact:1,chart:1,check:[1,2],chemic:1,children:1,choos:1,chosen:1,circumradiu:1,ckdtree:1,ckdtree_nod:1,ckdtreenod:1,classmethod:[1,2],cleanup:0,cleanup_build_fil:1,cleanup_core_fil:1,cleanup_result_fil:1,cleanup_tempfil:1,clear:1,close:1,closest:1,closest_vertex:1,cloud:1,cluster:1,code:4,coeffic:1,coeffici:1,collabor:1,colormap:1,column:1,com:[1,4],combin:1,compact:1,compact_nod:1,compar:[1,2],comparison:2,compil:[1,2],compile_prep:1,compos:1,comput:1,conain:1,conatain:1,concentr:1,condit:[1,2],config:1,connect:1,connolli:1,consid:1,consol:2,constant:[1,2],construct:1,constructor:2,constuct:1,contact:4,contain:[1,2],content:4,context:2,contian:1,contigu:1,continu:4,convers:2,convert:[1,2,3],convert_str:2,convet:2,coo_matrix:1,coodin:1,coordin:1,copi:[0,1,2,3],copy_data:1,copy_numb:1,copyright:[0,1,2,3],core:[0,4],corner:1,correl:1,correspond:[1,2],corrupt:1,cosmolog:1,could:[1,2],count:1,count_neighbor:1,coupl:1,cppconvert:2,cppexpressiontransform:2,cpu:1,cradded:1,creat:[1,2],create_2d_domain:1,create_3d_domain:1,csv:1,cumul:1,current:[1,2],custom:1,cylind:4,d:1,data:[1,2],data_funct:1,datafunct:0,datafunctionerror:1,dd:1,debug:[1,2],debug_level:[1,2],debugg:2,decreas:1,default_rng:1,defer:1,defin:[1,2],definit:1,degener:1,delete_all_paramet:1,delete_all_react:1,delete_all_speci:1,delete_paramet:1,delete_react:1,delete_speci:1,deltai:1,deltax:1,deltaz:1,demo:4,densiti:1,depend:1,dependency_graph:1,deprec:1,describ:[1,4],descript:1,detail:[0,1,2,3],determinist:[0,1,2,3,4],determinst:1,develop:[0,1,2,3,4],diamet:1,dict:[1,2],dictionari:[1,2],differ:[1,2],diffus:[0,1,2,3,4],diffusion_coeffici:1,dimens:1,dimension:1,direct:1,directli:[1,4],directori:1,discuss:1,displai:1,dissip:4,distanc:1,distance_between_2_vertic:1,distance_matrix:1,distance_upper_bound:1,distribut:[0,1,2,3],div:2,doe:1,doesn:1,dok_matrix:1,dolfin:1,dom:1,domain:[0,4],domainerror:1,domn:1,doubl:1,download:1,dr:4,draw:1,drawert:4,drawn:1,dtype:1,dual:1,duplic:2,durat:1,dure:[1,2],dynam:4,e:1,each:1,edg:1,edit:1,edu:1,effect:1,effici:1,either:1,element:1,ellipsi:2,els:1,embed:4,empti:2,end:1,engin:1,ensur:[1,2],enterpris:1,entir:1,entireti:1,ep:1,eq:2,equival:2,error:1,estim:1,et:1,euclidean:1,evalu:1,evaul:1,even:[0,1,2,3],evenli:1,exampl:1,exce:1,exceed:2,except:1,execut:[1,2],exist:1,expect:[1,2],expens:1,explor:1,export_stochss:1,export_to_csv:1,expres:2,express:[1,2],expressionconvert:2,expressionresult:2,extens:1,extra:1,f:1,f_durat:1,fail:2,fals:[1,2,3],far:2,faster:1,featur:4,fenic:1,field:1,fifo:1,figsiz:1,figur:1,figure_or_data:1,file:1,filenam:[1,3],fill:[1,4],fill_with_particl:1,find:[1,4],find_boundary_point:1,find_h:1,finit:1,first:[1,2],fit:[0,1,2,3],fix:1,flag:[1,2],flatten:1,flight:2,floordiv:2,flow:4,fluid:[1,4],folder_nam:1,follow:[1,2],forc:1,form:[1,4],format:1,forward:2,found:[1,2,4],foundat:[0,1,2,3],frac:1,fraction:1,frame:[1,4],free:[0,1,2,3],frequenc:1,from:[1,2],fromkei:1,fund:4,further:1,furthest:1,g:1,gain:1,galaxi:1,gdbgui:2,gener:[0,1,2,3,4],geometr:1,geometri:0,geometry_ivar:1,geometryal:1,geometryerror:1,geometryexterior:1,geometryinterior:1,get:[1,2],get_all_paramet:1,get_all_react:1,get_all_speci:1,get_arrai:1,get_array_nam:1,get_bounding_box:1,get_domain_s:1,get_num_point:1,get_num_react:1,get_num_speci:1,get_num_voxel:1,get_paramet:1,get_point:1,get_properti:1,get_react:1,get_speci:1,get_str:2,get_timespan:1,get_type_def:1,get_vol:1,getexpr_:2,getexpr_cpp:2,getexpr_python:2,github:[1,4],give:1,given:[1,2],gmsh:1,gnbu:1,gnu:[0,1,2,3,4],gprof:2,grai:1,grant:4,granular:2,graph:1,graph_obj:1,graph_refer:1,graviti:[1,4],greater:1,greatest:1,green:1,grei:1,gt:2,gte:2,guarante:1,h:1,ha:[1,4],hand:2,handl:2,happen:2,harvard:1,hasn:[1,2],have:[0,1,2,3],hawkin:1,height:1,held:1,help:1,here:[1,4],high:1,higher:1,highest:1,hit:1,hope:[0,1,2,3],how:[1,2,4],http:[0,1,2,3,4],hybrid:4,hyperrectangl:1,i:1,id:1,idea:1,ideal:2,idempot:1,identifi:[1,2],idx:1,ignor:1,ii:1,imag:1,image_height:1,image_width:1,immobil:1,implemenet:1,implement:[1,2],impli:[0,1,2,3],import_meshio_object:1,imposs:2,inch:1,incid:1,includ:1,included_types_list:1,incompat:1,increas:1,increment:2,indec:1,indeci:1,index:[1,4],indic:[1,2],indici:1,inf:1,infer:2,infin:1,infinit:1,info:1,inform:[1,2,4],inherit:2,init_cond:1,init_notebook_mod:1,initailcondit:1,initi:1,initialcondit:0,initialconditionerror:1,inlin:1,innernod:1,input:1,insert:1,insid:1,instanc:[1,2],instead:1,integ:1,intend:4,intern:1,internet:1,interv:1,invalid:[1,2],invalid_:2,invalid_nam:2,invalid_oper:2,io:1,iplot:1,ipython:1,is_valid:2,issuecom:1,item:[1,2],iter:1,its:1,j:1,jacob:4,jet:1,job:1,journal:1,jpeg:1,js:1,json:1,jupyt:1,k:1,kd:1,kd_tree1:1,kd_tree2:1,kd_tree:1,kdtree:1,kei:[1,2],keyerror:1,keyword:1,kill:1,kth:1,kwarg:1,l_i:1,lack:1,lagrangian:4,landi:1,larg:1,larger:1,last:[1,2],later:[1,2],latest:4,layout:1,leafnod:1,leafsiz:1,learn:1,leav:1,left:2,legaci:1,len:1,length:1,less:1,level:[1,2],li:4,librari:1,licens:[0,1,2,3,4],lifo:1,like:1,limit:1,linearli:1,link:1,link_text:1,list:[1,2],listofparamet:1,listofreact:1,listofspeci:1,listoftyp:1,load:1,local:1,locat:1,logic:2,longer:1,look:1,lookup:1,loos:1,lower:1,lowest:1,lshift:2,lt:2,lte:2,ly:1,m:1,mai:1,main:4,make:1,manag:2,maneewongvatana:1,manhattan:1,mani:[1,2],map:[1,2],map_oper:2,margin:1,mark:1,markers:1,mass:1,master:1,mathemat:[1,2],matmult:2,matplotlib:1,matrix:1,matter:1,max:1,max_dist:1,maximum:[1,2],mean:1,measur:1,median:1,meet:1,megabyt:1,merchant:[0,1,2,3],mesh:1,mesh_obj:1,meshio:1,messag:1,meta:1,method:[1,2],mgrid:1,midpoint:1,min:1,mine:1,minim:1,minimum:1,minkowski:1,miss:1,mod:2,model:[0,2,3,4],modelerror:1,modelname_mesh:1,modelname_species_:1,modifi:[0,1,2,3],modul:4,monthli:1,moor:1,more:[0,1,2,3,4],most:[1,2],mount:1,move:1,move_to_end:1,movement:1,msh:1,much:1,mult:2,multi:1,multipl:[1,2],must:1,n:1,n_i:1,naiv:1,name:[1,2],nameconst:2,namespac:[1,2],ndarrai:1,ndarri:1,nearbi:1,nearby_point:1,nearer:1,nearest:1,necessari:1,need:[1,4],neg:1,neighbor:1,neighbour:1,nih:4,node:[1,2],nodetransform:2,nodevisitor:2,non:1,none:[1,2,3],nonneg:1,norm:1,notat:1,note:[1,2],notebook:[1,4],noteq:2,notic:1,np:1,nsm:4,nu:1,num:2,num_step:1,number:[1,2],number_of_thread:[1,2],number_of_trajectori:[1,2],numpi:1,numpoint:1,nx:1,ny:1,nz:1,o:1,obj:1,object:[1,2],occur:1,od:1,off:1,offlin:1,og:1,omit:1,on_boundari:1,one:[1,2],onli:[1,2],onlin:1,open:1,oper:[1,2],operator_map:2,optim:1,option:1,orang:1,order:1,ordereddict:1,org:[0,1,2,3],origin:2,orrd:1,other:[1,2],otherwis:[1,2],our:[1,4],out:[1,4],outdat:1,output:[1,2],output_interv:1,output_typ:1,outsid:1,over:1,overflow:1,overrid:2,overridden:1,p0:1,p:1,p_ndx:1,packag:4,page:[1,4],pair:1,parallel:1,param:1,param_nam:1,param_obj:1,paramet:[0,2,3],parametererror:1,paramt:1,parent:1,pars:2,parse_comparison:2,parse_log:2,parse_oper:2,part:1,partcil:1,particip:1,particl:[1,4],particular:[0,1,2,3],pass:[1,2],path:[1,3],pdf:1,peculiar:1,per:1,perform:1,petzold:4,ph:1,physic:4,pii:4,pixel:1,place:1,placeinitialcondit:1,platform:1,pleas:[1,4],plot:1,plot_imag:1,plot_properti:1,plot_speci:1,plot_typ:1,plotli:1,plotly3:1,plt:1,pname:1,png:1,point:1,points1:1,points2:1,pop:[1,2],popitem:1,popul:1,posit:1,positon:1,pow:2,pre:2,predict:1,preper:1,present:[1,2],pressur:1,print:[1,2],problem:[0,1,2,3,4],process:[1,2],processor:1,produc:1,product:1,profil:[1,2],program:[0,1,2,3],project:4,propens:1,propensity_funct:1,properli:1,properti:[1,2],property_nam:1,proprti:1,protect:1,prove:4,provid:[1,2,4],prune:1,publish:[0,1,2,3],pubu:1,pubugn:1,pull:1,purd:1,purpl:1,purpos:[0,1,2,3],put:1,px:1,py:1,pypi:4,pyplot:1,python:[0,1,2,3,4],pythonconvert:2,pyurdm:4,quantifi:1,queri:1,query_ball_point:1,query_ball_tre:1,query_pair:1,question:4,quick:1,r:1,radii:1,radiu:1,rais:[1,2],random:[1,2],randomli:1,rang:1,rapidli:1,rate:1,rather:1,ratio:1,ravel:1,raw:1,rdme:4,rdpu:1,re:1,reac:1,reach:1,react_nam:1,react_obj:1,reactant:1,reaction:[0,2,3,4],reactionerror:1,read:1,read_fil:1,read_msh_fil:1,read_step:1,read_stochss_domain:1,read_stochss_subdomain_fil:1,read_xml_mesh:1,reader:1,real:1,receiv:[0,1,2,3],recent:[1,2],record:1,red:1,redistribut:[0,1,2,3],redshift:1,refer:1,referenc:1,regardless:1,region:1,regist:4,registr:4,reject:2,rememb:1,remov:1,replac:4,repres:1,represent:[1,2],request:1,requir:1,reserved_nam:1,resolut:1,restrict:1,restrict_to:1,result:[0,2],result_dir:1,resulterror:1,retreiv:1,retriev:1,return_length:1,return_plotly_figur:1,return_sort:1,return_stochss_model:[1,3],rho0:1,rho:1,right:[1,2],rname:1,rng:1,robust:1,royal:1,rshift:2,rtype:2,rule:1,run:[1,2],rvalu:2,s0021999118307101:4,s:[1,4],same:1,sampl:1,sanit:2,sanitized_data_function_nam:1,sanitized_parameter_nam:1,sanitized_species_nam:1,save:1,scalar:1,scale:1,scatter:1,scatterinitialcondit:1,schedul:1,scheme:1,scienc:[1,4],sciencedirect:4,scipi:1,sdm:1,sdpd:[1,4],search:[1,4],second:[1,2],section:1,see:[0,1,2,3,4],seed:[1,2],self:[1,2],sens:1,sep:1,separ:1,sequenc:1,seri:[1,2],set:[1,2],set_diffusion_coeffici:1,set_filenam:1,set_properti:1,set_timestep:1,setdefault:1,setter:1,shall:1,shallow:1,shape:1,sheth:1,should:[0,1,2,3],show:1,show_link:1,shrink:1,signal:1,significantli:1,similar:1,simlat:1,simplifi:1,simul:[0,1,2,3,4],simulationerror:[1,2],simulationtimeout:[1,2],sinc:1,singl:[1,2],size:1,skibba:1,sky:1,slide:1,slightli:1,small:1,smaller:1,smdl:1,smooth:4,sname:1,so:[1,2],societi:1,softwar:[0,1,2,3,4],solver:[0,1,4],some:1,sort:1,sound:1,sourc:[1,2,3,4],space:1,spactialpi:2,span:1,spars:1,sparse_distance_matrix:1,spatial:[0,1,2,3,4],spatialpy_build:1,spatialpy_model:[1,3],spatialpy_result:1,spatialpyerror:0,spatialpysolv:1,speci:0,special_charact:1,species_nam:1,species_object:1,specieserror:1,specif:1,specifi:1,speed:1,split:1,split_dim:1,squeez:1,src:1,ssa:1,sssa:4,start:[1,2,4],statement:2,statist:1,step:1,step_num:1,stochast:[0,1,2,3,4],stochss:[0,1,4],stochss_export:0,stoichiometr:1,stoichiometri:1,str:[1,2,3],string:[1,2],struct:2,structur:1,style:1,sub:2,subclass:[1,2],subdomain:1,submodul:0,subpackag:4,substanti:1,substitut:2,success:1,suggest:4,sum:1,suppli:1,support:[1,4],sure:1,survei:1,svg:1,symbol:2,syntax:2,syntaxerror:2,system:[1,4],szalai:1,t:[1,2],t_durat:1,t_ndx:1,t_ndx_list:1,t_val:1,tag:1,take:1,target:[1,2],tempfil:1,term:[0,1,2,3],tetrahedron:1,text:1,th:1,than:1,them:1,thi:[0,1,2,3,4],thin:1,those:1,though:1,thread:[1,2],time:[1,2],time_span:1,timeout:[1,2],timepoint:1,timespan:1,timestep:1,timestep_s:1,titl:1,toarrai:1,token:2,topolog:1,toroid:1,total:1,toth:1,trace:1,tracker:4,trajectori:[1,2],transit:1,travers:1,tree:[1,2],triangl:1,trigger:2,tspan:1,tupl:1,two:[1,2],txt:1,type:[1,2],type_id:1,u0:1,unabl:1,unaryop:2,under:[0,1,2,3,4],uniform:1,uniforminitialcondit:1,union:1,univers:[1,2],unless:1,unnecessari:1,unsqueez:1,unweight:1,up:1,updat:1,update_namespac:1,upload:1,upon:1,upper:1,urdm:1,us:[0,1,2,3,4],use_matplotlib:1,user:[1,4],usual:1,usub:2,util:1,v1:1,v:1,valid:[1,2],validationvisitor:2,valu:[1,2],valueerror:1,vari:1,varianc:1,veloc:1,verbos:2,version:[0,1,2,3,4],vertex:1,vertic:1,view:1,viewer:1,viscos:1,visit:2,visit_add:2,visit_and:2,visit_assign:2,visit_binop:2,visit_boolop:2,visit_byt:2,visit_cal:2,visit_compar:2,visit_const:2,visit_div:2,visit_ellipsi:2,visit_eq:2,visit_gt:2,visit_lt:2,visit_mult:2,visit_nam:2,visit_nameconst:2,visit_noteq:2,visit_num:2,visit_or:2,visit_pow:2,visit_str:2,visit_sub:2,visit_unaryop:2,visit_usub:2,visitor:2,vol:1,volum:1,voxel:1,vtk:1,vtkreader:0,vtkreadererror:1,vtkreaderioerror:1,wa:1,want:1,warranti:[0,1,2,3],we:[1,4],web:1,webp:1,weight:1,weir:4,well:2,were:[1,2],what:2,when:[1,2],where:1,whether:[1,2,3],which:[1,2],whose:1,wide:1,width:1,with_blacklist:2,with_namespac:2,within:1,without:[0,1,2,3],work:1,worker:1,wp:1,wrap:1,write:1,write_imag:1,written:1,www:[0,1,2,3,4],x1:1,x2:1,x:1,x_i:1,xaxi:1,xi:1,xk:1,xlim:1,xmax:1,xmin:1,xml:1,y:1,yaxi:1,yet:1,yi:4,ylgn:1,ylgnbu:1,ylim:1,ylorbr:1,ylorrd:1,ymax:1,ymin:1,you:[0,1,2,3,4],your:[1,4],z:1,zero:1,zeroth:1,zlim:1,zmax:1,zmin:1},titles:["spatialpy package","spatialpy.core package","spatialpy.solvers package","spatialpy.stochss package","Documentation for SpatialPy 1.0"],titleterms:{"0":4,"1":4,api:4,boundarycondit:1,build_express:2,cleanup:1,content:[0,1,2,3],copi:4,core:1,datafunct:1,document:4,domain:1,exampl:4,geometri:1,get:4,indic:4,initialcondit:1,issu:4,model:1,modul:[0,1,2,3],packag:[0,1,2,3],paramet:1,reaction:1,refer:4,report:4,result:1,solver:2,spatialpi:[0,1,2,3,4],spatialpyerror:1,speci:1,stochss:3,stochss_export:3,submodul:[1,2,3],subpackag:0,tabl:4,vtkreader:1}})