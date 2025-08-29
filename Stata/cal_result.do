clear
import excel data_od_vc_fdi_ind,firstrow
destring,replace
// 导入数据，并将字符转换为数字

forvalues x=2008/2018 {
	gen dummy`x'=0
	replace dummy`x'=1 if year==`x' 
}
//添加时间虚拟变量

gen TOTAL_E=solid+solid_cut
gen huji=pop_huji

drop if gdp==.
gen loggdp=log(gdp)
drop if labor==.
gen loglab=log(labor)
drop if capital==.
gen logcap=log(capital)
drop if TOTAL_E==.
gen logpoll=log(TOTAL_E)
drop if o_meanyear_==.
gen logomc=log(o_meanyear_)
drop if d_meanyear_==.
gen logdmc=log(d_meanyear_)
drop if FDI==.
gen logfdi=log(FDI)


reg loggdp loglab logcap logpoll
predict logA,residuals

//global unit=_b[_cons]
//global a_lab=_b[loglab]/(_b[loglab]+_b[logcap]+_b[logpoll])
//display $a_lab
//global a_cap=_b[logcap]/(_b[loglab]+_b[logcap]+_b[logpoll])
//display $a_cap
//global a_poll=1-$a_lab-$a_cap
//display $a_poll
//global r=0.20
//global gamma=1.3867
//global phi = 0.3510
//global psi = 1.5247
//计算产出函数里面三种要素的指数

//gen A = gdp/(labor^$a_lab)/(capital^$a_cap)/(TOTAL_E^$a_poll)
//gen A = exp(loggdp-$a_lab*loglab-$a_cap*logcap-$a_poll*logpoll)
//gen logA = log(A)
gen A = exp(logA)


gen logpollgdp = log(TOTAL_E/pop)-log(gdp/pop)
reg logpollgdp logomc logdmc logfdi dummy*
global m=exp(-_b[_cons])
display $m
predict reslogpopgdp, residuals
gen rho = exp(reslogpopgdp)/pm25

keep areacode year rho A pop huji city bl bc
reshape wide  rho* A* pop* huji* city* bl* bc*, i(areacode) j(year)
egen pop=rmean(pop*)
egen rho=rmean(rho*)
egen A=rmean(A*)
egen huji=rmean(huji*)
egen bl=rmean(bl*)
egen bc=rmean(bc*)
drop if pop==.
drop if A==.
drop if rho==.
drop if huji==.
gen lnrho=ln(rho)




keep areacode pop A rho huji bl bc
export excel data_ctrl, replace