# This file was generated, do not modify it. # hide
#hideall
@variables S_0 r d σ T
tgt_var=r
tgt_var2=σ
max_order=2
max_order2=1
# expr_el=generate_latex_expression_base(S_0,r,d,σ,T,tgt_var,max_order+max_order2)
expr_el=generate_latex_expression_base_mix(S_0,r,d,σ,T,tgt_var,max_order,tgt_var2,max_order2)
println(expr_el)