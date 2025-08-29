clear
import excel data_od_vc_fdi_ind,firstrow
save tempdata,replace
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

global unit=_b[_cons]
global r=0.20
global gamma=1.3867
global phi = 1/0.3510
global psi = 1.5247
//计算产出函数里面三种要素的指数


gen A = exp(logA)
drop if pop==.
gen logpop = log(pop)
gen logm=log((pop-huji)/huji)
gen logpm = log(pm25)
reg logpop logpm dummy*
reg logpop logA logm logomc logdmc logfdi dummy*
predict reslogpop, residuals
// 估计式：N-A

gen prelogpop = _b[logA]*logA+_b[logm]*logm
gen logpergdp = log(gdp/exp(logpop))
gen logperpoll = log(TOTAL_E/exp(prelogpop))
reg logperpoll logpergdp dummy*
//估计式：E-Y

gen prelogperpoll = _b[logpergdp]*logpergdp

predict reslogperpoll,residuals
gen respop_perpoll = reslogpop*reslogperpoll
reg logpop logA logm logomc logdmc logfdi reslogperpoll respop_perpoll dummy*
//分离特征


gen logperpollution=log(TOTAL_E/exp(logpop))
reg logperpollution logpergdp dummy*
test logpergdp-1=0
