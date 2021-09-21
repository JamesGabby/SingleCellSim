// Copyright (c) 2011-2015 by Thomas O'Hara, Yoram Rudy, 
//                            Washington University in St. Louis.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright 
// notice, this list of conditions and the following disclaimer in the 
// documentation and/or other materials provided with the distribution.
// 
// 3. Neither the names of the copyright holders nor the names of its
// contributors may be used to endorse or promote products derived from 
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF 
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND 
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// TypeScript translation of the C++ Implementation of the O'Hara-Rudy dynamic (ORd)
// model for the undiseased human ventricular action potential and calcium transient
//
// The ORd model is described in the article "Simulation of the Undiseased
// Human Cardiac Ventricular Action Potential: Model Formulation and
// Experimental Validation"
// by Thomas O'Hara, Laszlo Virag, Andras Varro, and Yoram Rudy
//
// The article and supplemental materails are freely available in the
// Open Access jounal PLoS Computational Biology
// Link to Article:
// http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1002061
// 
// Email: tom.ohara@gmail.com / rudy@wustl.edu
// Web: http://rudylab.wustl.edu
// 
revpots(); //compute reversal potentials //doubles = f64
RGC(); //compute rates, gates, and currents
stimulus(); //determine the value for the periodic stimulus
voltage(); //calculate the new membrane voltage
dVdt_APD(); //caluculate voltage derivative and APD90
FBC(); //calculate fluxes, buffers, and concentrations
var CL = 1000; //pacing cycle length
var ft = 1000 * CL; //final time
var skip = 10; //number of timesteps to skip in sampling of data in output file
var safetime = 25.0; //time from the beginning of each beat during which dt is fixed to small values
var beatssave = 2; //number of beats to save in the output
var amp = -80; //stimulus amplitude in uA/uF
var start = 0; //start time of the stimulus, relative to each beat		
var duration = 0.5; //duration of the stimulus in ms
var celltype = 0; //endo = 0, epi = 1, M = 2
//initial values for state variables, there are 41 of them
var v = -87.5;
var nai = 7;
var nass = nai;
var ki = 145;
var kss = ki;
var cai = 1.0e-4;
var cass = cai;
var cansr = 1.2;
var cajsr = cansr;
var m = 0;
var hf = 1;
var hs = 1;
var j = 1;
var hsp = 1;
var jp = 1;
var mL = 0;
var hL = 1;
var hLp = 1;
var a = 0;
var iF = 1;
var iS = 1;
var ap = 0;
var iFp = 1;
var iSp = 1;
var d = 0;
var ff = 1;
var fs = 1;
var fcaf = 1;
var fcas = 1;
var jca = 1;
var nca = 0;
var ffp = 1;
var fcafp = 1;
var xrf = 0;
var xrs = 0;
var xs1 = 0;
var xs2 = 0;
var xk1 = 1;
var Jrelnp = 0;
var Jrelp = 0;
var CaMKt = 0;
//constants
var nao = 140.0; //extracellular sodium in mM
var cao = 1.8; //extracellular calcium in mM
var ko = 5.4; //extracellular potassium in mM
//buffer paramaters
var BSRmax = 0.047;
var KmBSR = 0.00087;
var BSLmax = 1.124;
var KmBSL = 0.0087;
var cmdnmax = 0.05;
var kmcmdn = 0.00238;
var trpnmax = 0.07;
var kmtrpn = 0.0005;
var csqnmax = 10.0;
var kmcsqn = 0.8;
//CaMK paramaters
var aCaMK = 0.05;
var bCaMK = 0.00068;
var CaMKo = 0.05;
var KmCaM = 0.0015;
var KmCaMK = 0.15;
//physical constants
var R = 8314.0;
var T = 310.0;
var F = 96485.0;
//cell geometry
var L = 0.01;
var rad = 0.0011;
var vcell = 1000 * 3.14 * rad * rad * L;
var Ageo = 2 * 3.14 * rad * rad + 2 * 3.14 * rad * L;
var Acap = 2 * Ageo;
var vmyo = 0.68 * vcell;
var vmito = 0.26 * vcell;
var vsr = 0.06 * vcell;
var vnsr = 0.0552 * vcell;
var vjsr = 0.0048 * vcell;
var vss = 0.02 * vcell;
//introduce varaibles for reversal potentials, currents, fluxes, and CaMK
var ENa, EK, EKs;
var INa, INaL, Ito, ICaL, ICaNa, ICaK, IKr, IKs, IK1, INaCa_i, INaCa_ss, INaCa, INaK, IKb, INab, IpCa, ICab, Ist;
var Jrel, Jup, Jtr, Jdiff, JdiffNa, JdiffK, Jleak;
var CaMKa, CaMKb;
//introduce APD, timing, and counting parameters
var APD_flag = 0;
var APD;
var t_vdot_max;
var vrest;
var vo = v;
var dt = 0.005;
var t0 = 0;
var t = 0;
var dto;
var vdot_old;
var vdot = 0;
var vdot_max;
var p = 1;
var n = 0;
var count = 1;
//value holders for state varaibles in the case that the increase in dt was too aggressive, so a smaller one can be taken 
var nai0, nass0, ki0, kss0, cai0, cass0, cansr0, cajsr0, m0, hf0, hs0, jO, hsp0, jp0, mL0, hL0, hLp0, a0, iF0, iS0, ap0, iFp0, iSp0, d0, ff0, fs0, fcaf0, fcas0, jca0, nca0, ffp0, fcafp0, xrf0, xrs0, xs10, xs20, xk10, Jrelnp0, Jrelp0, CaMKt0;
var main = function () {
    // ouput data to file "output.txt" !! Must be used in Internet Explorer !!
    // let fso = new ActiveXObject("Scripting.FileSystemObject");
    // let address = "C:\\Users\\gabbs\\ord_model\\output.txt";
    // fso.CreateTextFile(address);
    // let file = fso.GetFile(address);
    // let streamWrite = file.OpenAsTextStream(2); //open the output file for appending new data
    // streamWrite.WriteLine("t\t v\t nai\t nass\t ki\t kss\t cai\t cass\t cansr\t cajsr\t Jrel\t CaMKt\t Jup\t Jtr\t Jdiff\t JdiffNa\t JdiffK\t Jleak\t INa\t INaL\t Ito\t ICaL\t ICaNa\t ICaK\t IKr\t IKs\t IK1\t INaCa_i\t INaCa_ss\t INaCa\t INaK\t IKb\t INab\t IpCa\t ICab\t Ist\t dt\t APD\t");
    console.log("Simulation running . . . ");
    while (t <= ft) {
        //rules for dynamic dt choice, and model letegration, comment to use fixed time steps
        if ((t >= (start + n * CL - 2) && t < (start + duration + n * CL)) ||
            (n >= 1 && t < (start + duration + (n - 1) * CL + safetime)) ||
            (APD_flag == 1 && v < 0.7 * vrest)) {
            dt = 0.005;
            t = t + dt;
            revpots();
            RGC();
            stimulus();
            vo = v;
            voltage();
            dVdt_APD();
            FBC();
        }
        else if (Math.abs(v - vo) < 0.2) {
            dt = Math.abs(0.8 / vdot);
            if (dt > 1.0) {
                dt = 1.0;
            }
            t = t + dt;
            revpots();
            RGC();
            stimulus();
            vo = v;
            voltage();
            dVdt_APD();
            FBC();
        }
        else if (Math.abs(v - vo) > 0.8) {
            nai0 = nai;
            nass0 = nass;
            ki0 = ki;
            cai0 = cai;
            cass0 = cass;
            cansr0 = cansr;
            cajsr0 = cajsr;
            m0 = m;
            hf0 = hf;
            hs0 = hs;
            jO = j;
            hsp0 = hsp;
            jp0 = jp;
            mL0 = mL;
            hL0 = hL;
            hLp0 = hLp;
            a0 = a;
            iF0 = iF;
            iS0 = iS;
            ap0 = ap;
            iFp0 = iFp;
            iSp0 = iSp;
            d0 = d;
            ff0 = ff;
            fs0 = fs;
            fcaf0 = fcaf;
            fcas0 = fcas;
            jca0 = jca;
            nca0 = nca;
            ffp0 = ffp;
            fcafp = fcafp;
            xrf0 = xrf;
            xrs0 = xrs;
            xs10 = xs1;
            xs20 = xs2;
            xk10 = xk1;
            Jrelnp0 = Jrelnp;
            Jrelp0 = Jrelp;
            CaMKt0 = CaMKt;
            t0 = t;
            dto = dt;
            dt = Math.abs(0.2 / vdot);
            t = t + dt;
            revpots();
            RGC();
            stimulus();
            vo = v;
            voltage();
            dVdt_APD();
            FBC();
            while (Math.abs(v - vo) > 0.8) {
                v = vo;
                nai = nai0;
                nass = nass0;
                ki = ki0;
                cai = cai0;
                cass = cass0;
                cansr = cansr0;
                cajsr = cajsr0;
                m = m0;
                hf = hf0;
                hs = hs0;
                j = jO;
                hsp = hsp0;
                jp = jp0;
                mL = mL0;
                hL = hL0;
                hLp = hLp0;
                a = a0;
                iF = iF0;
                iS = iS0;
                ap = ap0;
                iFp = iFp0;
                iSp = iSp0;
                d = d0;
                ff = ff0;
                fs = fs0;
                fcaf = fcaf0;
                fcas = fcas0;
                jca = jca0;
                nca = nca0;
                ffp = ffp0;
                fcafp = fcafp0;
                xrf = xrf0;
                xrs = xrs0;
                xs1 = xs10;
                xs2 = xs20;
                xk1 = xk10;
                Jrelnp = Jrelnp0;
                Jrelp = Jrelp0;
                CaMKt = CaMKt0;
                if (p == 1) {
                    dt = dto - 0.01;
                    p = 0;
                }
                else {
                    dt = dt - 0.01;
                }
                if (dt <= 0) {
                    dt = 1e-6;
                }
                t = t0 + dt;
                revpots();
                RGC();
                stimulus();
                voltage();
                dVdt_APD();
                FBC();
            }
            p = 1;
        }
        else {
            t = t + dt;
            revpots();
            RGC();
            stimulus();
            vo = v;
            voltage();
            dVdt_APD();
            FBC();
        }
        //uncomment below, and comment above to use a fixed dt
        /*t=t+dt; //fixed time step
        revpots();
        RGC();
        stimulus();
        vo=v;
        voltage();
        dVdt_APD();
        FBC();*/
        if (count % 500000 == 0) {
            console.log((t / ft * 100).toFixed(2) + "% complete"); //output runtime progress to the screen
        }
        //save results to output file when the sampling leterval and time are correct
        if (count % skip == 0 && t >= ft - beatssave * CL) {
            // 	fprintf(output,"%-018e\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\t%-07g\n",
            // 	t-(ft-beatssave*CL),v,nai,nass,ki,kss,cai,cass,cansr,cajsr,Jrel,CaMKt,Jup,Jtr,Jdiff,JdiffNa,JdiffK,Jleak,INa,INaL,Ito,ICaL,ICaNa,ICaK,IKr,IKs,IK1,INaCa_i,INaCa_ss,INaCa,INaK,IKb,INab,IpCa,ICab,Ist,dt,APD);
            //streamWrite.WriteLine(`${t-(ft-beatssave*CL)}\t ${v}\t ${nai}\t ${nass}\t ${ki}\t ${kss}\t ${cai}\t ${cass}\t ${cansr}\t ${cajsr}\t ${Jrel}\t ${CaMKt}\t ${Jup}\t ${Jtr}\t ${Jdiff}\t ${JdiffNa}\t ${JdiffK}\t ${Jleak}\t ${INa}\t ${INaL}\t ${Ito}\t ${ICaL}\t ${ICaNa}\t ${ICaK}\t ${IKr}\t ${IKs}\t ${IK1}\t ${INaCa_i}\t ${INaCa_ss}\t ${INaCa}\t ${INaK}\t ${IKb}\t ${INab}\t ${IpCa}\t ${ICab}\t ${Ist}\t ${dt}\t ${APD}\t`);
        }
        count++; //increase the loop counter
    }
    // fclose(output);//close the output file
    //streamWrite.close();
    console.log("Simulation complete");
    return 0;
};
function revpots() {
    ENa = (R * T / F) * Math.log(nao / nai);
    EK = (R * T / F) * Math.log(ko / ki);
    EKs = (R * T / F) * Math.log((ko + 0.01833 * nao) / (ki + 0.01833 * nai));
}
function RGC() {
    CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass);
    CaMKa = CaMKb + CaMKt;
    var vffrt = v * F * F / (R * T);
    var vfrt = v * F / (R * T);
    var mss = 1.0 / (1.0 + Math.exp((-(v + 39.57)) / 9.871));
    var tm = 1.0 / (6.765 * Math.exp((v + 11.64) / 34.77) + 8.552 * Math.exp(-(v + 77.42) / 5.955));
    m = mss - (mss - m) * Math.exp(-dt / tm);
    var hss = 1.0 / (1 + Math.exp((v + 82.90) / 6.086));
    var thf = 1.0 / (1.432e-5 * Math.exp(-(v + 1.196) / 6.285) + 6.149 * Math.exp((v + 0.5096) / 20.27));
    var ths = 1.0 / (0.009794 * Math.exp(-(v + 17.95) / 28.05) + 0.3343 * Math.exp((v + 5.730) / 56.66));
    var Ahf = 0.99;
    var Ahs = 1.0 - Ahf;
    hf = hss - (hss - hf) * Math.exp(-dt / thf);
    hs = hss - (hss - hs) * Math.exp(-dt / ths);
    var h = Ahf * hf + Ahs * hs;
    var jss = hss;
    var tj = 2.038 + 1.0 / (0.02136 * Math.exp(-(v + 100.6) / 8.281) + 0.3052 * Math.exp((v + 0.9941) / 38.45));
    j = jss - (jss - j) * Math.exp(-dt / tj);
    var hssp = 1.0 / (1 + Math.exp((v + 89.1) / 6.086));
    var thsp = 3.0 * ths;
    hsp = hssp - (hssp - hsp) * Math.exp(-dt / thsp);
    var hp = Ahf * hf + Ahs * hsp;
    var tjp = 1.46 * tj;
    jp = jss - (jss - jp) * Math.exp(-dt / tjp);
    var GNa = 75;
    var fINap = (1.0 / (1.0 + KmCaMK / CaMKa));
    INa = GNa * (v - ENa) * m * m * m * ((1.0 - fINap) * h * j + fINap * hp * jp);
    var mLss = 1.0 / (1.0 + Math.exp((-(v + 42.85)) / 5.264));
    var tmL = tm;
    mL = mLss - (mLss - mL) * Math.exp(-dt / tmL);
    var hLss = 1.0 / (1.0 + Math.exp((v + 87.61) / 7.488));
    var thL = 200.0;
    hL = hLss - (hLss - hL) * Math.exp(-dt / thL);
    var hLssp = 1.0 / (1.0 + Math.exp((v + 93.81) / 7.488));
    var thLp = 3.0 * thL;
    hLp = hLssp - (hLssp - hLp) * Math.exp(-dt / thLp);
    var GNaL = 0.0075;
    if (celltype == 1) {
        GNaL *= 0.6;
    }
    var fINaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
    INaL = GNaL * (v - ENa) * mL * ((1.0 - fINaLp) * hL + fINaLp * hLp);
    var ass = 1.0 / (1.0 + Math.exp((-(v - 14.34)) / 14.82));
    var ta = 1.0515 / (1.0 / (1.2089 * (1.0 + Math.exp(-(v - 18.4099) / 29.3814))) + 3.5 / (1.0 + Math.exp((v + 100.0) / 29.3814)));
    a = ass - (ass - a) * Math.exp(-dt / ta);
    var iss = 1.0 / (1.0 + Math.exp((v + 43.94) / 5.711));
    var delta_epi;
    if (celltype == 1) {
        delta_epi = 1.0 - (0.95 / (1.0 + Math.exp((v + 70.0) / 5.0)));
    }
    else {
        delta_epi = 1.0;
    }
    var tiF = 4.562 + 1 / (0.3933 * Math.exp((-(v + 100.0)) / 100.0) + 0.08004 * Math.exp((v + 50.0) / 16.59));
    var tiS = 23.62 + 1 / (0.001416 * Math.exp((-(v + 96.52)) / 59.05) + 1.780e-8 * Math.exp((v + 114.1) / 8.079));
    tiF *= delta_epi;
    tiS *= delta_epi;
    var AiF = 1.0 / (1.0 + Math.exp((v - 213.6) / 151.2));
    var AiS = 1.0 - AiF;
    iF = iss - (iss - iF) * Math.exp(-dt / tiF);
    iS = iss - (iss - iS) * Math.exp(-dt / tiS);
    var i = AiF * iF + AiS * iS;
    var assp = 1.0 / (1.0 + Math.exp((-(v - 24.34)) / 14.82));
    ap = assp - (assp - ap) * Math.exp(-dt / ta);
    var dti_develop = 1.354 + 1.0e-4 / (Math.exp((v - 167.4) / 15.89) + Math.exp(-(v - 12.23) / 0.2154));
    var dti_recover = 1.0 - 0.5 / (1.0 + Math.exp((v + 70.0) / 20.0));
    var tiFp = dti_develop * dti_recover * tiF;
    var tiSp = dti_develop * dti_recover * tiS;
    iFp = iss - (iss - iFp) * Math.exp(-dt / tiFp);
    iSp = iss - (iss - iSp) * Math.exp(-dt / tiSp);
    var ip = AiF * iFp + AiS * iSp;
    var Gto = 0.02;
    if (celltype == 1) {
        Gto *= 4.0;
    }
    if (celltype == 2) {
        Gto *= 4.0;
    }
    var fItop = (1.0 / (1.0 + KmCaMK / CaMKa));
    Ito = Gto * (v - EK) * ((1.0 - fItop) * a * i + fItop * ap * ip);
    var dss = 1.0 / (1.0 + Math.exp((-(v + 3.940)) / 4.230));
    var td = 0.6 + 1.0 / (Math.exp(-0.05 * (v + 6.0)) + Math.exp(0.09 * (v + 14.0)));
    d = dss - (dss - d) * Math.exp(-dt / td);
    var fss = 1.0 / (1.0 + Math.exp((v + 19.58) / 3.696));
    var tff = 7.0 + 1.0 / (0.0045 * Math.exp(-(v + 20.0) / 10.0) + 0.0045 * Math.exp((v + 20.0) / 10.0));
    var tfs = 1000.0 + 1.0 / (0.000035 * Math.exp(-(v + 5.0) / 4.0) + 0.000035 * Math.exp((v + 5.0) / 6.0));
    var Aff = 0.6;
    var Afs = 1.0 - Aff;
    ff = fss - (fss - ff) * Math.exp(-dt / tff);
    fs = fss - (fss - fs) * Math.exp(-dt / tfs);
    var f = Aff * ff + Afs * fs;
    var fcass = fss;
    var tfcaf = 7.0 + 1.0 / (0.04 * Math.exp(-(v - 4.0) / 7.0) + 0.04 * Math.exp((v - 4.0) / 7.0));
    var tfcas = 100.0 + 1.0 / (0.00012 * Math.exp(-v / 3.0) + 0.00012 * Math.exp(v / 7.0));
    var Afcaf = 0.3 + 0.6 / (1.0 + Math.exp((v - 10.0) / 10.0));
    var Afcas = 1.0 - Afcaf;
    fcaf = fcass - (fcass - fcaf) * Math.exp(-dt / tfcaf);
    fcas = fcass - (fcass - fcas) * Math.exp(-dt / tfcas);
    var fca = Afcaf * fcaf + Afcas * fcas;
    var tjca = 75.0;
    jca = fcass - (fcass - jca) * Math.exp(-dt / tjca);
    var tffp = 2.5 * tff;
    ffp = fss - (fss - ffp) * Math.exp(-dt / tffp);
    var fp = Aff * ffp + Afs * fs;
    var tfcafp = 2.5 * tfcaf;
    fcafp = fcass - (fcass - fcafp) * Math.exp(-dt / tfcafp);
    var fcap = Afcaf * fcafp + Afcas * fcas;
    var Kmn = 0.002;
    var k2n = 1000.0;
    var km2n = jca * 1.0;
    var anca = 1.0 / (k2n / km2n + Math.pow(1.0 + Kmn / cass, 4.0));
    nca = anca * k2n / km2n - (anca * k2n / km2n - nca) * Math.exp(-km2n * dt);
    var PhiCaL = 4.0 * vffrt * (cass * Math.exp(2.0 * vfrt) - 0.341 * cao) / (Math.exp(2.0 * vfrt) - 1.0);
    var PhiCaNa = 1.0 * vffrt * (0.75 * nass * Math.exp(1.0 * vfrt) - 0.75 * nao) / (Math.exp(1.0 * vfrt) - 1.0);
    var PhiCaK = 1.0 * vffrt * (0.75 * kss * Math.exp(1.0 * vfrt) - 0.75 * ko) / (Math.exp(1.0 * vfrt) - 1.0);
    var zca = 2.0;
    var PCa = 0.0001;
    if (celltype == 1) {
        PCa *= 1.2;
    }
    if (celltype == 2) {
        PCa *= 2.5;
    }
    var PCap = 1.1 * PCa;
    var PCaNa = 0.00125 * PCa;
    var PCaK = 3.574e-4 * PCa;
    var PCaNap = 0.00125 * PCap;
    var PCaKp = 3.574e-4 * PCap;
    var fICaLp = (1.0 / (1.0 + KmCaMK / CaMKa));
    ICaL = (1.0 - fICaLp) * PCa * PhiCaL * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCap * PhiCaL * d * (fp * (1.0 - nca) + jca * fcap * nca);
    ICaNa = (1.0 - fICaLp) * PCaNa * PhiCaNa * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaNap * PhiCaNa * d * (fp * (1.0 - nca) + jca * fcap * nca);
    ICaK = (1.0 - fICaLp) * PCaK * PhiCaK * d * (f * (1.0 - nca) + jca * fca * nca) + fICaLp * PCaKp * PhiCaK * d * (fp * (1.0 - nca) + jca * fcap * nca);
    var xrss = 1.0 / (1.0 + Math.exp((-(v + 8.337)) / 6.789));
    var txrf = 12.98 + 1.0 / (0.3652 * Math.exp((v - 31.66) / 3.869) + 4.123e-5 * Math.exp((-(v - 47.78)) / 20.38));
    var txrs = 1.865 + 1.0 / (0.06629 * Math.exp((v - 34.70) / 7.355) + 1.128e-5 * Math.exp((-(v - 29.74)) / 25.94));
    var Axrf = 1.0 / (1.0 + Math.exp((v + 54.81) / 38.21));
    var Axrs = 1.0 - Axrf;
    xrf = xrss - (xrss - xrf) * Math.exp(-dt / txrf);
    xrs = xrss - (xrss - xrs) * Math.exp(-dt / txrs);
    var xr = Axrf * xrf + Axrs * xrs;
    var rkr = 1.0 / (1.0 + Math.exp((v + 55.0) / 75.0)) * 1.0 / (1.0 + Math.exp((v - 10.0) / 30.0));
    var GKr = 0.046;
    if (celltype == 1) {
        GKr *= 1.3;
    }
    if (celltype == 2) {
        GKr *= 0.8;
    }
    IKr = GKr * Math.sqrt(ko / 5.4) * xr * rkr * (v - EK);
    var xs1ss = 1.0 / (1.0 + Math.exp((-(v + 11.60)) / 8.932));
    var txs1 = 817.3 + 1.0 / (2.326e-4 * Math.exp((v + 48.28) / 17.80) + 0.001292 * Math.exp((-(v + 210.0)) / 230.0));
    xs1 = xs1ss - (xs1ss - xs1) * Math.exp(-dt / txs1);
    var xs2ss = xs1ss;
    var txs2 = 1.0 / (0.01 * Math.exp((v - 50.0) / 20.0) + 0.0193 * Math.exp((-(v + 66.54)) / 31.0));
    xs2 = xs2ss - (xs2ss - xs2) * Math.exp(-dt / txs2);
    var KsCa = 1.0 + 0.6 / (1.0 + Math.pow(3.8e-5 / cai, 1.4));
    var GKs = 0.0034;
    if (celltype == 1) {
        GKs *= 1.4;
    }
    IKs = GKs * KsCa * xs1 * xs2 * (v - EKs);
    var xk1ss = 1.0 / (1.0 + Math.exp(-(v + 2.5538 * ko + 144.59) / (1.5692 * ko + 3.8115)));
    var txk1 = 122.2 / (Math.exp((-(v + 127.2)) / 20.36) + Math.exp((v + 236.8) / 69.33));
    xk1 = xk1ss - (xk1ss - xk1) * Math.exp(-dt / txk1);
    var rk1 = 1.0 / (1.0 + Math.exp((v + 105.8 - 2.6 * ko) / 9.493));
    var GK1 = 0.1908;
    if (celltype == 1) {
        GK1 *= 1.2;
    }
    if (celltype == 2) {
        GK1 *= 1.3;
    }
    IK1 = GK1 * Math.sqrt(ko) * rk1 * xk1 * (v - EK);
    var kna1 = 15.0;
    var kna2 = 5.0;
    var kna3 = 88.12;
    var kasymm = 12.5;
    var wna = 6.0e4;
    var wca = 6.0e4;
    var wnaca = 5.0e3;
    var kcaon = 1.5e6;
    var kcaoff = 5.0e3;
    var qna = 0.5224;
    var qca = 0.1670;
    var hca = Math.exp((qca * v * F) / (R * T));
    var hna = Math.exp((qna * v * F) / (R * T));
    var h1 = 1 + nai / kna3 * (1 + hna);
    var h2 = (nai * hna) / (kna3 * h1);
    var h3 = 1.0 / h1;
    var h4 = 1.0 + nai / kna1 * (1 + nai / kna2);
    var h5 = nai * nai / (h4 * kna1 * kna2);
    var h6 = 1.0 / h4;
    var h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna);
    var h8 = nao / (kna3 * hna * h7);
    var h9 = 1.0 / h7;
    var h10 = kasymm + 1.0 + nao / kna1 * (1.0 + nao / kna2);
    var h11 = nao * nao / (h10 * kna1 * kna2);
    var h12 = 1.0 / h10;
    var k1 = h12 * cao * kcaon;
    var k2 = kcaoff;
    var k3p = h9 * wca;
    var k3pp = h8 * wnaca;
    var k3 = k3p + k3pp;
    var k4p = h3 * wca / hca;
    var k4pp = h2 * wnaca;
    var k4 = k4p + k4pp;
    var k5 = kcaoff;
    var k6 = h6 * cai * kcaon;
    var k7 = h5 * h2 * wna;
    var k8 = h8 * h11 * wna;
    var x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    var x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    var x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    var x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    var E1 = x1 / (x1 + x2 + x3 + x4);
    var E2 = x2 / (x1 + x2 + x3 + x4);
    var E3 = x3 / (x1 + x2 + x3 + x4);
    var E4 = x4 / (x1 + x2 + x3 + x4);
    var KmCaAct = 150.0e-6;
    var allo = 1.0 / (1.0 + Math.pow(KmCaAct / cai, 2.0));
    var zna = 1.0;
    var JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    var JncxCa = E2 * k2 - E1 * k1;
    var Gncx = 0.0008;
    if (celltype == 1) {
        Gncx *= 1.1;
    }
    if (celltype == 2) {
        Gncx *= 1.4;
    }
    INaCa_i = 0.8 * Gncx * allo * (zna * JncxNa + zca * JncxCa);
    h1 = 1 + nass / kna3 * (1 + hna);
    h2 = (nass * hna) / (kna3 * h1);
    h3 = 1.0 / h1;
    h4 = 1.0 + nass / kna1 * (1 + nass / kna2);
    h5 = nass * nass / (h4 * kna1 * kna2);
    h6 = 1.0 / h4;
    h7 = 1.0 + nao / kna3 * (1.0 + 1.0 / hna);
    h8 = nao / (kna3 * hna * h7);
    h9 = 1.0 / h7;
    h10 = kasymm + 1.0 + nao / kna1 * (1 + nao / kna2);
    h11 = nao * nao / (h10 * kna1 * kna2);
    h12 = 1.0 / h10;
    k1 = h12 * cao * kcaon;
    k2 = kcaoff;
    k3p = h9 * wca;
    k3pp = h8 * wnaca;
    k3 = k3p + k3pp;
    k4p = h3 * wca / hca;
    k4pp = h2 * wnaca;
    k4 = k4p + k4pp;
    k5 = kcaoff;
    k6 = h6 * cass * kcaon;
    k7 = h5 * h2 * wna;
    k8 = h8 * h11 * wna;
    x1 = k2 * k4 * (k7 + k6) + k5 * k7 * (k2 + k3);
    x2 = k1 * k7 * (k4 + k5) + k4 * k6 * (k1 + k8);
    x3 = k1 * k3 * (k7 + k6) + k8 * k6 * (k2 + k3);
    x4 = k2 * k8 * (k4 + k5) + k3 * k5 * (k1 + k8);
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    KmCaAct = 150.0e-6;
    allo = 1.0 / (1.0 + Math.pow(KmCaAct / cass, 2.0));
    JncxNa = 3.0 * (E4 * k7 - E1 * k8) + E3 * k4pp - E2 * k3pp;
    JncxCa = E2 * k2 - E1 * k1;
    INaCa_ss = 0.2 * Gncx * allo * (zna * JncxNa + zca * JncxCa);
    INaCa = INaCa_i + INaCa_ss;
    var k1p = 949.5;
    var k1m = 182.4;
    var k2p = 687.2;
    var k2m = 39.4;
    k3p = 1899.0;
    var k3m = 79300.0;
    k4p = 639.0;
    var k4m = 40.0;
    var Knai0 = 9.073;
    var Knao0 = 27.78;
    var delta = -0.1550;
    var Knai = Knai0 * Math.exp((delta * v * F) / (3.0 * R * T));
    var Knao = Knao0 * Math.exp(((1.0 - delta) * v * F) / (3.0 * R * T));
    var Kki = 0.5;
    var Kko = 0.3582;
    var MgADP = 0.05;
    var MgATP = 9.8;
    var Kmgatp = 1.698e-7;
    var H = 1.0e-7;
    var eP = 4.2;
    var Khp = 1.698e-7;
    var Knap = 224.0;
    var Kxkur = 292.0;
    var P = eP / (1.0 + H / Khp + nai / Knap + ki / Kxkur);
    var a1 = (k1p * Math.pow(nai / Knai, 3.0)) / (Math.pow(1.0 + nai / Knai, 3.0) + Math.pow(1.0 + ki / Kki, 2.0) - 1.0);
    var b1 = k1m * MgADP;
    var a2 = k2p;
    var b2 = (k2m * Math.pow(nao / Knao, 3.0)) / (Math.pow(1.0 + nao / Knao, 3.0) + Math.pow(1.0 + ko / Kko, 2.0) - 1.0);
    var a3 = (k3p * Math.pow(ko / Kko, 2.0)) / (Math.pow(1.0 + nao / Knao, 3.0) + Math.pow(1.0 + ko / Kko, 2.0) - 1.0);
    var b3 = (k3m * P * H) / (1.0 + MgATP / Kmgatp);
    var a4 = (k4p * MgATP / Kmgatp) / (1.0 + MgATP / Kmgatp);
    var b4 = (k4m * Math.pow(ki / Kki, 2.0)) / (Math.pow(1.0 + nai / Knai, 3.0) + Math.pow(1.0 + ki / Kki, 2.0) - 1.0);
    x1 = a4 * a1 * a2 + b2 * b4 * b3 + a2 * b4 * b3 + b3 * a1 * a2;
    x2 = b2 * b1 * b4 + a1 * a2 * a3 + a3 * b1 * b4 + a2 * a3 * b4;
    x3 = a2 * a3 * a4 + b3 * b2 * b1 + b2 * b1 * a4 + a3 * a4 * b1;
    x4 = b4 * b3 * b2 + a3 * a4 * a1 + b2 * a4 * a1 + b3 * b2 * a1;
    E1 = x1 / (x1 + x2 + x3 + x4);
    E2 = x2 / (x1 + x2 + x3 + x4);
    E3 = x3 / (x1 + x2 + x3 + x4);
    E4 = x4 / (x1 + x2 + x3 + x4);
    var zk = 1.0;
    var JnakNa = 3.0 * (E1 * a3 - E2 * b3);
    var JnakK = 2.0 * (E4 * b1 - E3 * a1);
    var Pnak = 30;
    if (celltype == 1) {
        Pnak *= 0.9;
    }
    if (celltype == 2) {
        Pnak *= 0.7;
    }
    INaK = Pnak * (zna * JnakNa + zk * JnakK);
    var xkb = 1.0 / (1.0 + Math.exp(-(v - 14.48) / 18.34));
    var GKb = 0.003;
    if (celltype == 1) {
        GKb *= 0.6;
    }
    IKb = GKb * xkb * (v - EK);
    var PNab = 3.75e-10;
    INab = PNab * vffrt * (nai * Math.exp(vfrt) - nao) / (Math.exp(vfrt) - 1.0);
    var PCab = 2.5e-8;
    ICab = PCab * 4.0 * vffrt * (cai * Math.exp(2.0 * vfrt) - 0.341 * cao) / (Math.exp(2.0 * vfrt) - 1.0);
    var GpCa = 0.0005;
    IpCa = GpCa * cai / (0.0005 + cai);
}
function FBC() {
    var CaMKb = CaMKo * (1.0 - CaMKt) / (1.0 + KmCaM / cass);
    CaMKa = CaMKb + CaMKt;
    CaMKt += dt * (aCaMK * CaMKb * (CaMKb + CaMKt) - bCaMK * CaMKt);
    JdiffNa = (nass - nai) / 2.0;
    JdiffK = (kss - ki) / 2.0;
    Jdiff = (cass - cai) / 0.2;
    var bt = 4.75;
    var a_rel = 0.5 * bt;
    var Jrel_inf = a_rel * (-ICaL) / (1.0 + Math.pow(1.5 / cajsr, 8.0));
    if (celltype == 2) {
        Jrel_inf *= 1.7;
    }
    var tau_rel = bt / (1.0 + 0.0123 / cajsr);
    if (tau_rel < 0.005) {
        tau_rel = 0.005;
    }
    Jrelnp = Jrel_inf - (Jrel_inf - Jrelnp) * Math.exp(-dt / tau_rel);
    var btp = 1.25 * bt;
    var a_relp = 0.5 * btp;
    var Jrel_infp = a_relp * (-ICaL) / (1.0 + Math.pow(1.5 / cajsr, 8.0));
    if (celltype == 2) {
        Jrel_infp *= 1.7;
    }
    var tau_relp = btp / (1.0 + 0.0123 / cajsr);
    if (tau_relp < 0.005) {
        tau_relp = 0.005;
    }
    Jrelp = Jrel_infp - (Jrel_infp - Jrelp) * Math.exp(-dt / tau_relp);
    var fJrelp = (1.0 / (1.0 + KmCaMK / CaMKa));
    Jrel = (1.0 - fJrelp) * Jrelnp + fJrelp * Jrelp;
    var Jupnp = 0.004375 * cai / (cai + 0.00092);
    var Jupp = 2.75 * 0.004375 * cai / (cai + 0.00092 - 0.00017);
    if (celltype == 1) {
        Jupnp *= 1.3;
        Jupp *= 1.3;
    }
    var fJupp = (1.0 / (1.0 + KmCaMK / CaMKa));
    Jleak = 0.0039375 * cansr / 15.0;
    Jup = (1.0 - fJupp) * Jupnp + fJupp * Jupp - Jleak;
    Jtr = (cansr - cajsr) / 100.0;
    nai += dt * (-(INa + INaL + 3.0 * INaCa_i + 3.0 * INaK + INab) * Acap / (F * vmyo) + JdiffNa * vss / vmyo);
    nass += dt * (-(ICaNa + 3.0 * INaCa_ss) * Acap / (F * vss) - JdiffNa);
    ki += dt * (-(Ito + IKr + IKs + IK1 + IKb + Ist - 2.0 * INaK) * Acap / (F * vmyo) + JdiffK * vss / vmyo);
    kss += dt * (-(ICaK) * Acap / (F * vss) - JdiffK);
    var Bcai;
    if (celltype == 1) {
        Bcai = 1.0 / (1.0 + 1.3 * cmdnmax * kmcmdn / Math.pow(kmcmdn + cai, 2.0) + trpnmax * kmtrpn / Math.pow(kmtrpn + cai, 2.0));
    }
    else {
        Bcai = 1.0 / (1.0 + cmdnmax * kmcmdn / Math.pow(kmcmdn + cai, 2.0) + trpnmax * kmtrpn / Math.pow(kmtrpn + cai, 2.0));
    }
    cai += dt * (Bcai * (-(IpCa + ICab - 2.0 * INaCa_i) * Acap / (2.0 * F * vmyo) - Jup * vnsr / vmyo + Jdiff * vss / vmyo));
    var Bcass = 1.0 / (1.0 + BSRmax * KmBSR / Math.pow(KmBSR + cass, 2.0) + BSLmax * KmBSL / Math.pow(KmBSL + cass, 2.0));
    cass += dt * (Bcass * (-(ICaL - 2.0 * INaCa_ss) * Acap / (2.0 * F * vss) + Jrel * vjsr / vss - Jdiff));
    cansr += dt * (Jup - Jtr * vjsr / vnsr);
    var Bcajsr = 1.0 / (1.0 + csqnmax * kmcsqn / Math.pow(kmcsqn + cajsr, 2.0));
    cajsr += dt * (Bcajsr * (Jtr - Jrel));
}
function voltage() {
    v += -dt * (INa + INaL + Ito + ICaL + ICaNa + ICaK + IKr + IKs + IK1 + INaCa + INaK + INab + IKb + IpCa + ICab + Ist);
}
function stimulus() {
    if ((t > (start + n * CL) && t < (start + duration + n * CL - dt))) {
        if (Ist == 0) {
            vrest = v;
        }
        Ist = amp;
    }
    else if (t > (start + duration + n * CL - dt)) {
        Ist = 0.0;
        n = n + 1;
    }
}
function dVdt_APD() {
    vdot_old = vdot;
    vdot = (v - vo) / dt;
    if (APD_flag == 0 && v > -40 && vdot < vdot_old) {
        vdot_max = vdot_old;
        t_vdot_max = t - dt;
        APD_flag = 1;
    }
    if (APD_flag == 1 && v < 0.9 * vrest) {
        APD = t - t_vdot_max;
        APD_flag = 0;
    }
}
