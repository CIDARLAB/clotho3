angular.module("clotho.dna",[]),angular.module("clotho.dna").service("DNA",[function(){var a={},b={},c={},d={},e={},f={},g={},h={};a.letters=/[A-Z]/gi,a.nucleotide=/[gautcGAUTC]/g,a.dna=/[gatcGATC]/g,a.rna=/[gaucGAUC]/g,a.nucleotide_degnerate=/[gatucryswkmbdhvnxGATUCRYSWKMBDHVNX\^_]/g,a.protein=/[ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g,a.protein_degenerate=/[ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g,a.strip={},a.strip.dna=/[^gatcGATC]/g,a.strip.nucleotide_degenerate=/[^gatucryswkmbdhvnxGATUCRYSWKMBDHVNX]/g,a.strip.rna=/[^gaucGAUC]/g,a.strip.protein=/[^ACDEFGHIKLMNPQRSTVWYZacdefghiklmnpqrstvwyz\*]/g,a.strip.protein_degenerate=/[^ABCDEFGHIKLMNPQRSTVWYXZabcdefghiklmnpqrstvwyxz\*]/g,a.strip.alignment=/[^\.\-]/g,a.strip.non_letters=/[^A-Z]/gi,a.strip.whitespace=/\s/g,b.dna="acgt".split(""),b.rna="acug".split(""),b.nucleotide="acgtu".split(""),b.nucleotide_degenerate="acgturyswkmbdhvnx".split(""),b.protein="ACDEFGHIKLMNPQRSTVWYZ*".split(""),b.protein_degnerate="ABCDEFGHIKLMNPQRSTVWXYZ*".split(""),b.all="ABCDEFGHIKLMNPQRSTUVWXYZ*".split(""),c.nucleotide_degenerate={A:"A",B:"[CGTU]",C:"C",D:"[AGTU]",G:"G",H:"[ACTU]",K:"[GTU]",M:"[AC]",N:"[ACGTU]",R:"[AG]",S:"[CG]",T:"T",U:"U",V:"[ACG]",W:"[ATU]",Y:"[CTU]",".":"[ACGTU]","-":"[ACGTU]"},c.nucleotide_undegenerate={A:"A",C:"C",G:"G",T:"T",U:"U","[ACGTU]":"N","[ACG]":"V","[ACTU]":"H","[AC]":"M","[AGTU]":"D","[AG]":"R","[ATU]":"W","[CGTU]":"B","[CG]":"S","[CTU]":"Y","[GTU]":"K"},c.dna_degenerate={A:"A",B:"[CGT]",C:"C",D:"[AGT]",G:"G",H:"[ACT]",K:"[GT]",M:"[AC]",N:"[ACGT]",R:"[AG]",S:"[CG]",T:"T",V:"[ACG]",W:"[AT]",Y:"[CT]",".":"[ACGT]","-":"[ACGT]"},c.rna_degenerate={A:"A",B:"[CGU]",C:"C",D:"[AGU]",G:"G",H:"[ACU]",K:"[GU]",M:"[AC]",N:"[ACGU]",R:"[AG]",S:"[CG]",U:"U",V:"[ACG]",W:"[AU]",Y:"[CU]",".":"[ACGU]","-":"[ACGU]"},c.amino_one_to_three={A:"Ala",B:"Asx",C:"Cys",D:"Asp",E:"Glu",F:"Phe",G:"Gly",H:"His",I:"Ile",K:"Lys",L:"Leu",M:"Met",N:"Asp",P:"Pro",Q:"Gln",R:"Arg",S:"Ser",T:"Thr",V:"Val",W:"Trp",X:"X",Y:"Tyr",Z:"Glx","*":"*"},c.amino_one_to_full={A:"Alanine",B:"[Aspartic Acid or Asparagine]",C:"Cysteine",D:"Aspartic Acid",E:"Glutamic Acid",F:"Phenylalanine",G:"Glycine",H:"Histidine",I:"Isoleucine",K:"Lysine",L:"Leucine",M:"Methionine",N:"Asparagine",P:"Proline",Q:"Glutamine",R:"Arginine",S:"Serine",T:"Threonine",V:"Valine",W:"Tryptophan",X:"Any",Y:"Tyrosine",Z:"[Glutamine or Glutamic Acid]","*":"Stop"},d.dna={a:"t",c:"g",g:"c",t:"a",A:"T",C:"G",G:"C",T:"A",r:"y",y:"r",R:"Y",Y:"R",k:"m",m:"k",K:"M",M:"K",b:"v",v:"b",B:"V",V:"B",d:"h",h:"d",D:"H",H:"D",n:"n",N:"N","^":"_",_:"^"},d.rna={a:"u",c:"g",g:"c",u:"a",A:"U",C:"G",G:"C",U:"A",r:"y",y:"r",R:"Y",Y:"R",k:"m",m:"k",K:"M",M:"K",b:"v",v:"b",B:"V",V:"B",d:"h",h:"d",D:"H",H:"D",n:"n",N:"N","^":"_",_:"^"},d.transcribe={a:"u",c:"g",g:"c",t:"a",A:"U",C:"G",G:"C",T:"A"},d.reverse_transcribe={a:"t",c:"g",g:"c",u:"a",A:"T",C:"G",G:"C",U:"A"},d.translate={uuu:"F",ucu:"S",uau:"Y",ugu:"C",uuc:"F",ucc:"S",uac:"Y",ugc:"C",uua:"L",uca:"S",uaa:"*",uga:"*",uug:"L",ucg:"S",uag:"*",ugg:"W",cuu:"L",ccu:"P",cau:"H",cgu:"R",cuc:"L",ccc:"P",cac:"H",cgc:"R",cua:"L",cca:"P",caa:"Q",cga:"R",cug:"L",ccg:"P",cag:"Q",cgg:"R",auu:"I",acu:"T",aau:"N",agu:"S",auc:"I",acc:"T",aac:"N",agc:"S",aua:"I",aca:"T",aaa:"K",aga:"R",aug:"M",acg:"T",aag:"K",agg:"R",guu:"V",gcu:"A",gau:"D",ggu:"G",guc:"V",gcc:"A",gac:"D",ggc:"G",gua:"V",gca:"A",gaa:"E",gga:"G",gug:"V",gcg:"A",gag:"E",ggg:"G"},d.reverse_translate={A:["gcu","gcc","gca","gcg"],C:["ugu","ugc"],D:["gau","gac"],E:["gaa","gag"],F:["uuu","uuc"],G:["ggu","ggc","gga","ggg"],H:["cau","cac"],I:["auu","auc","aua"],K:["aaa","aag"],L:["uua","uug","cuu","cuc","cua","cug"],M:["aug"],N:["aau","aac"],P:["ccu","ccc","cca","ccg"],Q:["caa","cag"],R:["cgu","cgc","cga","cgg","aga","agg"],S:["ucu","ucc","uca","ucg","agu","agc"],T:["acu","acc","aca","acg"],V:["guu","guc","gua","gug"],W:["ugg"],Y:["uau","uac"],"*":["uaa","uag","uga"]},d.reverse_translate_regexp={A:/gc[acgturyswkmbdhvn]/,C:/[tu]g[ctuy]/,D:/ga[tcuy]/,E:/ga[agr]/,F:/[tu][tu][tcuy]/,G:/gg[acgturyswkmbdhvn]/,H:/ca[tcuy]/,I:/a[tu][atcuwmhy]/,K:/aa[agr]/,L:/c[tu][acgturyswkmbdhvn]|[tu][tu][agr]|[ctuy][tu][agr]/,M:/a[tu]g/,N:/aa[tucy]/,P:/cc[acgturyswkmbdhvn]/,Q:/ca[agr]/,R:/cg[acgturyswkmbdhvn]|ag[agr]|[cam]g[agr]/,S:/[tu]c[acgturyswkmbdhvn]|ag[ct]/,T:/ac[acgturyswkmbdhvn]/,V:/g[tu][acgturyswkmbdhvn]/,W:/[tu]gg/,Y:/[tu]a[ctuy]/,"*":/[tu]a[agr]|[tu]ga|[tu][agr]a/},e.standard={},f.ecoli={},g.ds_terminal={g:-2.8,a:4.1,t:4.1,c:-2.8},g.dh_terminal={g:.1,a:2.3,t:2.3,c:.1},g.ds={gg:-19.9,ag:-21,tg:-22.7,cg:-27.2,ga:-22.2,aa:-22.2,ta:-21.3,ca:-22.7,gt:-22.4,at:-20.4,tt:-22.2,ct:-21,gc:-27.2,ac:-22.4,tc:-22.2,cc:-19.9},g.dh={gg:-8,ag:-7.8,tg:-8.5,cg:-10.6,ga:-8.2,aa:-7.9,ta:-7.2,ca:-8.5,gt:-8.4,at:-7.2,tt:-7.9,ct:-7.8,gc:-10.6,ac:-8.4,tc:-8.2,cc:-8},h.dna={a:313.2,t:304.2,c:289.2,g:329.2,u:290.2},h.rna={a:329.2,t:320.2,c:305.2,g:345.2,u:306.2};var i=function(b){return b.replace(a.letters,"")},j=function(b,c){var d=c?a.strip.dna:a.strip.nucleotide_degenerate;return b.replace(d,"")},k=function(a,b){return-1==a.search(b)},l=function(b,d){return d=d?d:c.nucleotide_degenerate,b.replace(a.nucleotide_degenerate,function(a){return d[a.toUpperCase()]})},m=function(a){return a.split("").reverse().join("")},n=function(b){return b.replace(a.nucleotide_degnerate,function(a){return d.dna[a]})},o=function(a){return a.replace(/[uU]/gi,function(a){return{u:"t",U:"T"}[a]})},p=function(a){return a.replace(/[tT]/gi,function(a){return{t:"u",T:"U"}[a]})},q=function(a){var b=[];return _.forEach(a.split(""),function(a){b.push(c.amino_one_to_three[a])}),b.join(" ")},r=function(a){var b=[];return _.forEach(a.split(""),function(a){b.push(c.amino_one_to_full[a])}),b.join(" ")},s=function(a){return m(n(a))},t=function(a,c){if(c=c||b.dna,"string"==typeof a&&(a=a.length),0==a)return"";for(var d=[],e=0,f="",g=0;a>g;g++)e=Math.floor(Math.random()*c.length),f=c[e],d.push(f);return d.join("")},u=function(a){return _.shuffle(a.split("")).join("")},v=function(a){return a.match(/.{1,3}/g).join(" ")},w=function(b){return b.replace(a.dna,function(a){return d.transcribe[a]})},x=function(b){return b.replace(a.rna,function(a){return d.reverse_transcribe[a]})},y=function(b,c){k(b,a.strip.dna)&&(b=p(b));var e=b.length;if(3>e)return"";b=b.toLowerCase();var f=c?c:B(b),g="";b=b.substring(f),e%3!=0&&(b=b.substring(0,3*Math.floor(e/3)));for(var h=0;h<b.length;h+=3)g+=d.translate[b.substr(h,3)];return g},z=function(a,b){if(a){for(var c=!b||_.isUndefined(b)||null==b?!0:!1,e="",f=a.split(""),g=0;g<f.length;g++){var h=d.reverse_translate[f[g]];e+=c?h[Math.floor(Math.random()*h.length)]:h[0]}return e}},A=function(a){var b={},c=new RegExp("((atg)(.+?)(ta[agr]|tga|t[agr]a))","gi");if(!c.test(a))return!1;for(var d=0;3>d;d++){b[d]={},b[d].sequence=a,b[d].frags=[];for(var e;null!=(e=c.exec(a));)b[d].frags.push({match:e[0],index:e.index,length:e[0].length});b[d].frags.sort(function(a,b){return b.length-a.length}),b[d].longest=b[d].frags[0].length}return b},B=function(a){var b=0,c=0,d=A(a);if(!d)return 0;for(var e=0;e<d.length;e++)c<d[e].longest&&(c=d[e].longest,b=e);return b},C=function(a,b,c){var d=c?l(b):b,e=a.match(new RegExp("(?=("+d+"))","gi"));return e?e.length:0},D=function(a,c){c="undefined"!=typeof c?c:b.all;var d={};return _.forEach(c,function(b){d[b]=C(a,b)}),d},E=function(a,c,d){d="undefined"!=typeof d?d:b.dna,c="undefined"!=typeof c?c:0;var e={};return _.forEach(d,function(b){_.forEach(d,function(d){var f=b+d,g=C(a,f);g>=c&&(e[f]=g)})}),e},F=function(a){var b=a.length,c=a.replace(/[^GCS]/gi,"").length;return c/b},G=function(a){{var c=0,d=0,e=0;D(a,b.nucleotide),E(a,b.nucleotide),F(a)}return[c,d,e]},H=function(a){if(a.length<14){var c=D(a,b.dna);return 2*(c.a+c.t)+4*(c.c+c.g)}return 64.9+41*(F(a)*a.length-16.4)/a.length},I=function(a,c){c=c?c:.05;var d=D(a,b.dna);return a.length<14?2*(d.a+d.t)+4*(d.c+d.g)-16.6*(Math.log(.05)/Math.log(10))+16.6*(Math.log(c)/Math.log(10)):a.length<50?100.5+41*(d.g+d.c)/(d.a+d.t+d.g+d.c)-820/(d.a+d.t+d.g+d.c)+16.6*(Math.log(c)/Math.log(10)):81.5+41*(d.g+d.c)/(d.a+d.t+d.g+d.c)-500/(d.a+d.t+d.g+d.c)+16.6*(Math.log(c)/Math.log(10))},J=function(b,c){if(!(b.length<1)){if(a.strip.dna.exec(b))return console.log("contains non-dna letters... aborting melting temp calc"),void 0;b=b.toLowerCase(),c="undefined"!=typeof c?c:{},c.dna="undefined"!=typeof c.dna?c.dna:2e-7,c.salt="undefined"!=typeof c.salt?c.salt:.05,c.mg="undefined"!=typeof c.mg?c.mg:.0015;var d=1.987,e=0,f=0,h=c.salt+140*c.mg;e+=.368*(b.length-1)*Math.log(h),e=e+g.ds_terminal[b.charAt(0)]+g.ds_terminal[b.charAt(b.length-1)],f=f+g.dh_terminal[b.charAt(0)]+g.dh_terminal[b.charAt(b.length-1)];var i=E(b);return _.forEach(i,function(a,b){e+=g.ds[b]*a,f+=g.dh[b]*a}),1e3*f/(e+d*Math.log(c.dna/2))-273.15}},K=function(c){if(c&&!(c.length<1)){var d=k(c,a.strip.dna)?"dna":"rna",e=b.nucleotide,f=0,g=D(c,e);return _.forEach(e,function(a){f+=h[d][a]*g[a]}),f+=159}},L=function(a,b){return new Array(+b+1).join(a)},M=function(a,c){var d=!0;c="undefined"!=typeof c?c:4,b="undefined"!=typeof b?b:b.nucleotide;var e=E(a,c);for(var f in e){if(!e.hasOwnProperty(f))break;new RegExp(L(f,c),"ig").test(a)&&(d=!1)}return d},N=function(a,c){var d=!0;c="undefined"!=typeof c?c:4,b=b.nucleotide;for(var e=0;e<b.length;e++)C(a,L(b[e],c))&&(d=!1);return d},O=function(a){for(var b=0,c=0,d=0;d<a.length;d++){var e=J(a[d]);b=b>e?b:e,c=e>c?c:e}return 6>=b-c};return{regexps:a,monomers:b,maps:c,complements:d,geneticCodes:e,frequencies:f,weights:h,lettersOnly:i,dnaOnly:j,verify:k,undegenerize:l,reverse:m,complement:n,rna_to_dna:o,dna_to_rna:p,amino_one_to_three:q,amino_one_to_full:r,revcomp:s,randomSequence:t,shuffleSequence:u,codonSpace:v,transcribe:w,reverseTranscribe:x,translate:y,reverseTranslate:z,calcORFs:A,probableORF:B,occuranceCount:C,monomer_count:D,neighbor_count:E,GC_content:F,gibbs:G,melting_temp_basic:H,melting_temp_saltAdjusted:I,melting_temp:J,molecular_weight:K,createRun:L,verifyRepeats:M,verifyRuns:N,verifyMeltingTemps:O}}]),angular.module("clotho.dna").directive("dnaShuffle",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.shuffleSequence)}}}]),angular.module("clotho.dna").directive("dnaRandom",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.randomSequence)}}}]),angular.module("clotho.dna").directive("dnaReverse",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.reverse)}}}]),angular.module("clotho.dna").directive("dnaComplement",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.complement)}}}]),angular.module("clotho.dna").directive("dnaRevcomp",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.revcomp)}}}]),angular.module("clotho.dna").directive("dnaTranscribe",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){var f=function(b){return console.log(b),a.transcribe(b)};e.$formatters.push(f)}}}]),angular.module("clotho.dna").directive("dnaReverseTranscribe",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){var f=function(b){return a.reverseTranscribe(b)};e.$formatters.push(f)}}}]),angular.module("clotho.dna").directive("dnaTranslate",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){var f=function(b){return a.translate(b)};e.$formatters.push(f)}}}]),angular.module("clotho.dna").directive("dnaReverseTranslate",["DNA","$parse",function(a,b){return{restrict:"A",require:"ngModel",link:function(c,d,e,f){var g=b(e.dnaReverseTranslate)(c)||!1,h=function(b){return a.reverseTranslate(b,g)};f.$formatters.push(h)}}}]),angular.module("clotho.dna").directive("dnaDnaOnly",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.dnaOnly)}}}]),angular.module("clotho.dna").directive("dnaRemoveNonLetters",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.lettersOnly)}}}]),angular.module("clotho.dna").directive("dnaAminoOneToThree",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.amino_one_to_three)}}}]),angular.module("clotho.dna").directive("dnaAminoOneToFull",["DNA",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$formatters.push(a.amino_one_to_full)}}}]),angular.module("clotho.dna").directive("dnaGcContent",["DNA","$filter",function(a,b){return{restrict:"A",require:"ngModel",link:function(c,d,e,f){var g=function(c){return b("number")(a.GC_content(c))};f.$formatters.push(g)}}}]),angular.module("clotho.dna").directive("dnaMeltingTemp",["DNA","$filter",function(a,b){return{restrict:"A",require:"ngModel",link:function(c,d,e,f){var g=function(c){return b("number")(a.melting_temp(c))};f.$formatters.push(g)}}}]),angular.module("clotho.dna").directive("dnaMolecularWeight",["DNA","$filter",function(a,b){return{restrict:"A",require:"ngModel",link:function(c,d,e,f){var g=function(c){return b("number")(a.molecular_weight(c))};f.$formatters.push(g)}}}]),angular.module("clotho.dna").service("Digest",["DNA",function(a){var b={BglII:{name:"BglII",match:"agatct",cut:"a^gatc_t",strand:"",methylation:!1,overhang:4,type:"II",subtype:"S",notes:{},buffer:"","star activity":!1,comment:"",rebase:"http://rebase.neb.com/rebase/enz/BglII.html",personal:{},citations:{},ordering:{}},BsaI:{name:"BsaI",match:"ggtctc",cut:"ggtctc (1/5)",strand:"",methylation:!1,overhang:3,type:"II",subtype:"P",notes:{},buffer:"","star activity":!1,comment:"",rebase:"http://rebase.neb.com/rebase/enz/BsaI.html",personal:{},citations:{},ordering:{}},BsmbI:{name:"BsmbI",match:"cgtctc",cut:"cgtctc (1/5)",strand:"",methylation:!0,overhang:4,type:"II",subtype:"S",notes:{},buffer:"","star activity":!1,comment:"",rebase:"http://rebase.neb.com/rebase/enz/BsmbI.html",personal:{},citations:{},ordering:{}},XhoI:{name:"XhoI",match:"ctcgag",cut:"c^tcga_g",strand:"",methylation:!0,overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!1,comment:"flanking C5 methylation can slow cleavage, ATCTCTCGAGTCTA is cut v. slowly",rebase:"http://rebase.neb.com/rebase/enz/XhoI.html",personal:{},citations:{},ordering:{}},BamHI:{name:"BamHI",match:"ggatcc",cut:"g^gatc_c",strand:"",methylation:{},overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!1,comment:"",rebase:"http://rebase.neb.com/rebase/enz/BamHI.html",personal:{},citations:{},ordering:{}},EcoRI:{name:"EcoRI",match:"gaattc",cut:"g^aatt_c",strand:"",methylation:{},overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!0,comment:"",rebase:"http://rebase.neb.com/rebase/enz/EcoRI.html",personal:{},citations:{},ordering:{}},XbaI:{name:"XbaI",match:"tctaga",cut:"t^ctag_a",strand:"",methylation:{},overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!0,comment:"",rebase:"http://rebase.neb.com/rebase/enz/XbaI.html",personal:{},citations:{},ordering:{}},SpeI:{name:"XbaI",match:"actagt",cut:"a^ctag_t",strand:"",methylation:!0,overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!0,comment:"",rebase:"http://rebase.neb.com/rebase/enz/SpeI.html",personal:{},citations:{},ordering:{}},PstI:{name:"PstI",match:"ctgcag",cut:"c_tgca^g",strand:"",methylation:{},overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!0,comment:"",rebase:"http://rebase.neb.com/rebase/enz/PstI.html",personal:{},citations:{},ordering:{}},HindIII:{name:"HindIII",match:"aagctt",cut:"a^agct_t",strand:"",methylation:{},overhang:4,type:"II",subtype:"P",notes:{},buffer:"","star activity":!0,comment:"",rebase:"http://rebase.neb.com/rebase/enz/HindIII.html",personal:{},citations:{},ordering:{}},AlwnI:{name:"AlwnI",match:"cagnnnctg",cut:"cag_nnn^ctg",strand:"",methylation:{},overhang:3,type:"II",subtype:"P",notes:{},buffer:"","star activity":!1,comment:"",rebase:"http://rebase.neb.com/rebase/enz/AlwnI.html",personal:{},citations:{},ordering:{}},AcuI:{name:"AcuI",match:"ctgaag",cut:"ctgaag (16/14)",strand:"",methylation:{},overhang:2,type:"II",subtype:"G, S, alpha",notes:{},buffer:"","star activity":!1,comment:"Methylated product: 6-methyladenosine (base undetermined)",rebase:"http://rebase.neb.com/rebase/enz/XbaI.html",personal:{},citations:{},ordering:{}},BseRI:{name:"BseRI",match:"gaggag",cut:"gaggag (10/8)",strand:"",methylation:{},overhang:2,type:"II",subtype:"G, S",notes:{},buffer:"","star activity":!1,comment:"Methylated Product: 6-methyladenosine (base undetermined)",rebase:"http://rebase.neb.com/rebase/enz/XbaI.html",personal:{},citations:{},ordering:{}}},c={};_.forEach(b,function(a,b){c[a.match]=b,c[a.cut]=b});var d={blunt:{mark:"|",type:"Blunt",description:'A cut that results in no "sticky" ends, or overhangs. '},main:{mark:"^",type:"Main Strand",description:"Denotes a cut on the 5' -> 3' (primary) strand, usually visualized as the \"top\" strand."},comp:{mark:"_",type:"Complementary Strand",description:"Denotes a cut on the 3' -> 5' (complementary) strand, usually visualized as the \"bottom\" strand."}},e={};e.localCuts=/\(((.*?)([\^|_])(.+?)([\^|_])(.*?))|((.*?)(\|)(.*?))\)/gi,e.nonLocalCuts=/\((\d+)\/(\d+)\)/gi,e.findCut=/(\|)|([\^_])(.+?)([_\^])/gi,e.findBlunt=/(\|)/gi,e.findOverhang=/([\^_])(.+?)([_\^])/gi;var f=function(a,b){return a+a.substr(0,b)},g=function(a){if(a.length>1){var b=a.shift();a[a.length-1]=a[a.length-1]+b}return a},h=function(b,c){var d=c?b+"|"+a.revcomp(b):b;return new RegExp(a.undegenerize(d),"ig")},i=function(a,b){return h(b,!0).test(f(a,b.length-1))},j=function(a){return c[a]},k=function(b,c,d,e,f){d=_.isUndefined(d)?3:d;var g=a.randomSequence(d)+(f?a.revcomp(c.match):c.match)+a.randomSequence(d);return e?g+b:b+g},l=function(b,c,d){var e=d?d:a.probableORF(b),f=m(b,c,!0),g=Object.keys(f);if(0==f.length)return b;for(var i=c instanceof RegExp?c:h(c,!0),j=0;j<g.length;j++){var k=g[j],l=f[k],n=k-k%3+e,o=Math.ceil((l.length+(k-n))%3);a:for(var p=0;o>=p;p++){var q=n+3*p,r=b.substr(q,3).toLowerCase(),s=a.complements.translate[r],t=a.complements.reverse_translate[s];if(!(t.length<=1)){for(var u=0;u<t.length;u++)if(t[u]!=r){var v=b.substr(0,q)+t[u]+b.substr(q+3);if(!i.test(v)){b=v;break a}}i.test()}}if(impossible)return!1}return b},m=function(a,b,c){for(var d,e={},g=b instanceof RegExp?b:h(b,c!==!1);null!=(d=g.exec(f(a,b.length-1)));)e[d.index]=b;return e},n=function(a,b,c){return Object.keys(m(a,b,c))},o=function(a){return a.replace(/[\|\^_]/gi,"")},p=function(a){return a.replace(/[\(\)]/gi,"")},q=function(a){return p(o(a))},r=function(b,c){return c?(c=_.isArray(c)?c:[c],_.each(c,function(c){var d=/\((\d+)\/(\d+)\)/gi.exec(c.cut);if(d){var e={};e.enz="("+a.undegenerize(c.match)+")",e.rev="("+a.undegenerize(a.revcomp(c.match))+")";var f=["^","_"],g=["_","^"],i=d[1]<d[2]?f:g;if(d[1]<0){e.gap1="(.{"+Math.abs(d[1]-d[2])+"})",e.gap2="(.{"+Math.abs(d[2])+"})";var j=new RegExp(e.gap2+e.gap1+e.enz,"ig");b=b.replace(j,function(a,b,c,d){return[i[0]+b+i[1]+c+"("+d+")"]});var k=new RegExp(e.enz+e.gap1+e.gap2,"ig");b=b.replace(k,function(a,b,c,d){return["("+b+")"+c+i[0]+d+i[1]]})}else{e.gap1="(.{"+d[1]+"})",e.gap2="(.{"+(d[2]-d[1])+"})";var j=new RegExp(e.enz+e.gap1+e.gap2,"ig");b=b.replace(j,function(a,b,c,d){return["("+b+")"+c+i[0]+d+i[1]]});var k=new RegExp(e.gap2+e.gap1+e.rev,"ig");b=b.replace(k,function(a,b,c,d){return[i[0]+b+i[1]+c+"("+d+")"]})}}else b=b.replace(h(c.match),"("+c.cut+")"),b=b.replace(h(a.revcomp(c.match)),"("+a.revcomp(c.cut)+")");var l=b.substring(b.length-c.match.length)+b.substr(c.match.length-1);h(c.match).exec(l)}),b):b},s=function(a,b){return o(r(a,b))},t=function(a,b){return p(r(a,b))},u=function(a,b){for(var c,d=e.findCut,f=[];null!=(c=d.exec(a));)c.match=c[0],c.isBlunt="|"==c.match,c.is3prime=c.isBlunt?null:"_"==c[2],c.overhang=c.isBlunt?"":c[3],c.length=c.match.length,c.terminal=!(0!=c.index&&c.index+(c.isBlunt?1:c.length)!=a.length),b?!c.terminal&&f.push(c):f.push(c);return f},v=function(a,b){var c=0,d="",e=[],f=[],h=u(a);for(var i in h)h.hasOwnProperty(i)&&e.push(i);e.sort();for(var j=0;j<e.length;j++){var k=h[e[j]],l=k.isBlunt?"":k.match,m=a.substring(c,k.index)+l;f.push(m),c=k.index,d=l}return f.push(a.substring(c)),b?g(f):f},w=function(a){return _.sortBy(a,function(a){return-a.length})},x=function(a,b){if(_.isString(a))return a;var c;return _.isUndefined(b)?(b=0,c=a.length-1):c=0,_.sortBy(a,function(a){return Math.abs(a.length-b)})[c]},y=function(a){return a.replace(/(_.+?\^)|(\^.+?_)|(\|)/gi,"")},z=function(a,b){var c=_.find(u(a),function(a){return console.log(a),!a.terminal});return b&&c.index<a.length-c.index-c.length?a.substring(c.index):a.substring(0,c.index+c.length)},A=function(a){if(a.indexOf("_")<0)return a;var b=/(.*?)(_.+?\^?.*)/gi,c=b.exec(a),d=(c[2],c[1]);return d},B=function(){},C=function(){},D=function(a,b,c,d){if(!b)return"no enzyme provided";d&&(a=d(a)),a=t(a,b);var e=v(a,c);return _.sortBy(e,function(a){return-a.length})};return{enzymes:b,cutMarks:d,extendSequence:f,createRegex:h,circularize:g,sitePresent:i,identifySite:j,findMatches:m,findIndices:n,markSites:r,markMatches:s,markCuts:t,removeCuts:o,removeMatches:p,removeMarks:q,findOverhangs:u,makeCuts:v,addRestrictionSite:k,swapoutSites:l,sortFragments:w,gelPurify:x,removeOverhangs:y,trimPastInternal:z,exonuclease35:A,exonuclease53:B,polymerase53:C,digest:D}}]),angular.module("clotho.dna").directive("digestMark",["Digest","$filter","$parse",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){var f;b.$watch(d.digestMark,function(a){f=a,e.$render()});var g=function(b){return a.markSites(b,f)},h=function(b){return a.removeMarks(b)};e.$parsers.unshift(h),e.$formatters.push(g)}}}]),angular.module("clotho.dna").directive("digestHighlight",["Digest","$parse","$compile","$filter",function(a,b,c,d){return{restrict:"A",require:"ngModel",link:function(b,e,f,g){function h(){var f=a.markSites(g.$modelValue,b.highlightEnz),h=/\((.+?)\)/gi,i=f.replace(h,"<digest-annotation>$1</digest-annotation>");i=d("DigestCuts")(i),e.html(c("<div>"+i+"</div>")(b))}b.$watch(f.digestHighlight,function(a){console.log(a),b.highlightEnz=a,h()}),b.$watch(function(){return g.$modelValue},function(){h()})}}}]),angular.module("clotho.dna").directive("digestAnnotation",["$tooltip",function(){return{restrict:"EA",replace:!1,transclude:!0,template:'<span tooltip="{{ highlightEnz.name }}" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',compile:function(){return{pre:function(){},post:function(a,b){b.css({backgroundColor:"#fcc"})}}}}}]),angular.module("clotho.dna").directive("digestCutTop",[function(){return{restrict:"EA",link:function(a,b){b.html("&#8595;"),b.css("color","#f00")}}}]),angular.module("clotho.dna").directive("digestCutBottom",[function(){return{restrict:"EA",link:function(a,b){b.html("&#8593;"),b.css("color","#00f")}}}]),angular.module("clotho.dna").filter("DigestCuts",function(){return function(a,b){var c=a.replace(/\^/gi,b?"":"<digest-cut-top></digest-cut-top>");return c=c.replace(/_/gi,b?"":"<digest-cut-bottom></digest-cut-bottom>")}}),angular.module("clotho.dna").service("PCR",["DNA","Digest",function(a,b){function c(a){a instanceof c&&(a=a.sequence),b.findOverhangs(a,!0).length,this.sequence=a,this.process()}c.prototype={process:function(a){a&&(this.sequence=a),this.ends=b.findOverhangs(this.sequence),this.endMatches=_.pluck(this.ends,"match")},reverse:function(){this.sequence=a.revcomp(this.sequence),this.process()},circularize:function(){2==this.endMatches.length&&this.endsMatch(this.endMatches[0],this.endMatches[1])&&(this.sequence=this.sequence.substring(0,this.sequence.length-this.endMatches[1].length),this.sequence=this.ends[0].overhang+this.sequence.substring(this.endMatches[0].length),this.process())},endsMatch:function(b,c){return b==c||b==a.revcomp(c)},endsMatchRevcomp:function(b,c){return b==a.revcomp(c)},findMatchingEnds:function(a){return _.filter(this.ends,function(b){return this.endsMatch(b.match,a)},this)},findMatchingEndIndices:function(a){var b=this.findMatchingEnds(a);switch(b.length){case 0:return!1;case 1:return[_.indexOf(this.ends,b[0])];default:return[0,1]}},canMatch:function(a){return a instanceof c?this.canMatchFrag(a):!!this.findMatchingEnds(a).length},canMatchFrag:function(a){return!!(_.find(a.endMatches,function(a){return this.canMatch(a)},this)||[]).length},matchMap:function(a){var b={};return _.each(this.endMatches,function(c,d){_.each(a.endMatches,function(a,e){this.endsMatch(c,a)&&(b[d]?console.log("pointing to two fragments"):b[d]=e)},this)},this),_.keys(b).length?b:void 0},matchMapArray:function(a){var b={};return _.each(a,function(a,c){this!==a&&(b[c]=this.matchMap(a))},this),b},findFirstMatch:function(a){return _.find(a,function(a){return this!==a&&this.canMatch(a)},this)},findFirstMatchIndex:function(a){return _.indexOf(a,this.findFirstMatch(a))},canMatchArray:function(a){return!!this.findFirstMatch(a)},joinFragment:function(a){_.isString(a)&&(a=new c(a));var d=this.matchMap(a);if(!_.keys(d).length)return!1;var e=_.pairs(d)[0],f=e[0],g=e[1],h=this.ends[f],i=a.ends[g];0==h.index?0==i.index&&(g=g||1==a.ends.length?0:1,a.reverse(),i=a.ends[g]):0!=i.index&&(g=g||1==a.ends.length?0:1,a.reverse(),i=a.ends[g]);var j=0==h.index?a:this,k=0==h.index?i:h,l=j===this?a:this,m=k===h?i:h,n=j.sequence.substring(0,k.index)+b.removeMarks(k.match)+l.sequence.substring(m.index+m.match.length);return this.process(n),!0},alignFragment:function(){}};var d=function(c,d){var e=b.findIndices(c,d,!1),f=b.findIndices(c,a.revcomp(d),!1);return{forward:e,reverse:f}},e=function(c,d){for(var e,f,g=null,h=8,i=d.length-h,j={};e=d.substring(i),f=b.createRegex(e);--i){if(console.log(i),j.forward=c.match(f)||[],j.reverse=a.revcomp(c).match(f)||[],console.log(j,j.forward.length,j.reverse.length,j.forward.length+j.reverse.length),!j.forward.length&&!j.reverse.length){console.log("no *exact* matches found for length "+i+" from 3 prime end"),g=!1;break}if(j.forward.length+j.reverse.length==1){console.log("one match"),g=j.forward.length?c.search(f):a.revcomp(c).search(f);break}console.log("else")}return g},f=function(a,b,c){return c=!!c||!0,c?e(a,b):e(a,b)},g=function(a,b,c){if(_.isEmpty(b)||_.isEmpty(c))return"a primer is not defined";var e=d(a,b),f=d(a,c);return e.forward.length+e.reverse.length<1?"primer1 : no matches":f.forward.length+f.reverse.length<1?"primer2 : no matches":e.forward.length+e.reverse.length>1?"primer1 : multiple matches":f.forward.length+f.reverse.length>1?"primer2 : multiple matches":e.forward.length&&f.forward.length||e.reverse.length&&f.reverse.length?"primers point same direction":!0},h=function(a,b){if(2!=b.length)return"Can only handle having two primers right now";var c=b[0],e=b[1],f=g(a,c,e);if(f!==!0)return f;var h=d(a,c),k=d(a,e),l=h.forward.length?+h.forward[0]:+h.reverse[0],m=k.forward.length?+k.forward[0]:+k.reverse[0];return h.forward.length?m>l?i(a,l,m+e.length):j(a,l,m+e.length):l>m?i(a,m,l+c.length):j(a,m,l+c.length)},i=function(a,b,c){return a.substring(b,c)},j=function(a,b,c){return a.substring(b)+a.substring(0,c)},k=function(b,c){if(2!=c.length)return"Can only handle having two primers right now";var e=c[0],f=c[1],h=g(b,e,f);if(h!==!0)return h;var i=d(b,e),j=d(b,f),k=i.forward.length?+i.forward[0]:+i.reverse[0],l=j.forward.length?+j.forward[0]:+j.reverse[0],m=a.createRun(" ",k);return m+=i.forward.length?c[0]:a.revcomp(c[0]),m+=a.createRun(" ",l-k-c[0].length),m+=j.forward.length?c[1]:a.revcomp(c[1]),m+=a.createRun(" ",b.length-m.length)},l=function(a,d){var e=[];return _.each(a,function(a){if(_.isArray(a)&&(a=a[0]),b.findOverhangs(a,!0).length)if(d){var f=b.trimPastInternal(a,!0);e.push(new c(f))}else _.each(b.makeCuts(a),function(a){e.push(new c(a))});else e.push(new c(a))}),e},m=function(a){a=l(a,!0);for(var b=0;b<a.length&&1!=a.length;b++){var c=a[b],d=c.findFirstMatchIndex(a);d>=0&&(c.joinFragment(a[d]),a.splice(d,1),b--)}return _.each(a,function(a){a.circularize()}),1==a.length?a[0].sequence:_.pluck(a,"sequence")};return{predict:h,anneal:f,findAnnealSimple:e,primerAlign:k,parseFragments:l,ligate:m}}]),angular.module("clotho.dna").directive("pcrPredict",["PCR","Digest","DNA",function(a){return{restrict:"A",require:"ngModel",scope:{backbone:"=ngModel",primers:"="},link:function(b,c,d,e){function f(){c.text(a.predict(b.backbone,b.primers))}e.$render=function(){f()},b.$watch(function(){return b.primers[0]+b.primers[1]},function(){f()})}}}]),angular.module("clotho.dna").directive("pcrAlign",["PCR","Digest","DNA","$filter",function(a,b,c,d){return{restrict:"A",require:"ngModel",scope:{backbone:"=ngModel",primers:"="},link:function(b,c,e,f){function g(){var e=a.primerAlign(b.backbone,b.primers),f=57;e=d("breakLines")(e,f,"*").split("*");for(var g=d("breakLines")(b.backbone,f,"*").split("*"),h="",i=0;i<g.length;i++)h+=e[i]+"\n"+g[i]+"\n";c.html(h)}f.$render=function(){g()},b.$watch("primers",g,!0)}}}]),angular.module("clotho.dna").directive("ligateAlign",["PCR","Digest","DNA","$compile","$filter",function(a,b,c,d,e){return{restrict:"A",require:"ngModel",scope:{fragments:"=ngModel"},link:function(b,c,f,g){function h(){var f=a.ligate(b.fragments,!0,!0);console.log(f),f=_.isArray(f)?"did not ligate to completion : "+JSON.stringify(f):e("DigestCuts")(f,!0),c.html(d("<span>"+f+"</span>")(b))}g.$render=function(){h()},b.$watch("fragments",function(){h()},!0)}}}]),angular.module("clotho.dna").directive("ligateFrag",function(){return{restrict:"EA",replace:!1,transclude:!0,template:'<span tooltip="Initial Fragment" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',link:function(a,b){b.css("color","#faa")}}}),angular.module("clotho.dna").directive("ligateStickymatch",function(){return{restrict:"EA",replace:!1,transclude:!0,template:'<span tooltip="Sticky-end Complementary Region" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',link:function(a,b){b.css("color","#6b6")}}}),angular.module("clotho.dna").directive("ligateBluntmatch",function(){return{restrict:"EA",replace:!1,transclude:!0,template:'<span tooltip="Note! Blunt ends will only yield this product 50% of the time (fragment direction is random)" tooltip-placement="mouse" tooltip-animation="false" tooltip-append-to-body="true" ng-transclude></span>',link:function(a,b){b.css("color","#6b6")}}}),angular.module("clotho.dna").service("Construction",["DNA","Digest","PCR","$parse","$q",function(a,b,c,d,e){var f={};f.DNA=a,f.Digest=b,f.PCR=c;var g={"PCR.predict":{reaction:"PCR.predict",input:["",[]],readable:"PCR"},"PCR.ligate":{reaction:"PCR.ligate",input:[[]],readable:"Ligation"},"Digest.digest":{reaction:"Digest.digest",input:["",[]],readable:"Restriction Digest"},"Digest.gelPurify":{reaction:"Digest.gelPurify",input:[],readable:"Gel Purify"}},h=function(a,b){var c=a.dictionary;if(angular.isEmpty(c))return b?[]:"";console.log(c),a.dictionaryObject={};var g={};console.log("PROCESSING CONSTRUCTION FILE",a),_.remove(c,function(a){return!!a.computed});var h=[];angular.forEach(c,function(a){if(angular.isString(a.value)&&angular.extend(a,{computed:!1}),a.Clotho){if(a.preprocess&&(!a.retrieved||a.process!=a.preprocess)){var b=d(a.preprocess)(f);angular.isUndefined(b),h.push=angular.extend(a,{processed:a.preprocess,retrieved:!0,computed:!1,value:b})}}else a.preprocess&&(console.log(a.value),a.value=angular.isObject(a.value)?a.value:JSON.parse(JSON.stringify(a.value)))});var i=e.defer();return e.all(h).then(function(b){console.log(b),console.log(c),_.each(c,function(a){g[a.key]=a}),console.log(g);var h=e.when();return angular.forEach(a.steps,function(a,b){h.then(function(){var h=e.defer(),i=[];angular.forEach(a.input,function(a){if(angular.isString(a))b=d(a)(g)||"",b=b.value||b,i.push(b);else{var b=angular.isArray(a)?[]:{};angular.forEach(a,function(a){var c=d(a)(g)||{};c=c.value||c,b.push(c)}),i.push(b)}});var j=d(a.reaction)(f).apply(null,i),k={key:a.output,computed:!0,stepNum:b,value:j};
return c.push(k),g[k.key]=k,h.resolve(),h.promise})}),h}).then(function(){a.dictionaryObject=g,i.resolve(b?c:c.final.value)}),i.promise};return{reactions:g,process:h}}]),angular.module("clotho.dna").filter("orderByProp",function(){return function(a){return _.sortBy(a,function(a){return a.prop})}}),angular.module("clotho.dna").filter("objAddKeys",function(){return function(a){return _.each(a,function(a,b){_.extend(a,{key:b})})}}),angular.module("clotho.dna").filter("stepCheck",function(){return function(a,b){return _.pick(a,function(a){return 0==a.computed||a.stepNum<b})}}),angular.module("clotho.dna").directive("constructionDictionaryView",[function(){return{restrict:"EA",require:"ngModel",scope:{dictionary:"=ngModel",editable:"=constructionEditable",uptostep:"=constructionUptostep",processto:"=constructionProcessto"},templateUrl:"views/_dna/construction/dictionary.html",link:function(a){a.uptofilter=function(b){return!b.computed||b.stepNum<a.uptostep&&b.stepNum<a.processto},a.addTerm=function(){console.log(a.dictionary),a.dictionary.push({key:"",value:""})}}}}]),angular.module("clotho.dna").directive("constructionField",["$compile","$filter",function(a,b){return{restrict:"EA",require:"ngModel",scope:{fieldName:"@constructionField",model:"=ngModel",dictionary:"=?constructionDictionary",options:"=constructionOptions",editable:"=constructionEditable",stepIndex:"=constructionStepindex",required:"@constructionRequired",type:"@constructionType"},compile:function(){return{pre:function(b,c){function d(a){switch(angular.lowercase(a)){case"value":return'<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-disabled="!editable" ng-required="required"/>';case"boolean":return'<button btn-checkbox class="btn" ng-class="{\'btn-success\' : !!model}"  ng-model="model" ng-disabled="!editable"><i ng-class="{\'icon-white icon-ok-circle\' : !!model, \'icon-ban-circle\' : !model}"></i></button>';case"array":return angular.isString(b.model)&&(b.model=[b.model]),'<input class="span12" type="text" placeholder="{{fieldName}}" ng-model="model" ng-list ng-change="checkInput()" ng-disabled="!editable" ng-required="required"/>';default:return'<select class="span12" placeholder="{{fieldName}}" ng-model="model" ng-options="key as key for (key,val) in dict" ng-disabled="!editable" ng-change="checkInput()" ng-required="required"></select>'}}var e=angular.element('<div class="form-group" ng-class="{\'has-error\' : hasInvalidItem || hasFutureItem}"><label class="control-label">{{ fieldName }}</label><constructionFieldInput></constructionFieldInput></div>');e.find("constructionFieldInput").html(d(b.type)),c.html(a(e)(b))},post:function(a){a.required="false"!=a.required?"true":"false",a.$watch("dictionary",function(){a.dict=b("stepCheck")(a.dictionary,a.stepIndex)}),a.$watch("dict",function(){a.checkInput()},!0),a.checkInput=function(){function b(b){a.dict[b]?a.dict[b].stepNum>=a.stepIndex&&(a.hasFutureItem=!0):a.hasInvalidItem=!0}a.model&&"boolean"!=a.type&&"value"!=a.type&&(a.hasInvalidItem=!1,a.hasFutureItem=!1,"array"==a.type?angular.forEach(a.model,function(a){b(a)}):b(a.model))}}}}}}]),angular.module("clotho.dna").directive("constructionStep",["Construction","$parse","$compile","$http","$templateCache","$filter","$timeout","$dialog",function(a,b,c,d,e,f,g,h){return{restrict:"EA",require:"ngModel",scope:{step:"=ngModel",dictionary:"=constructionDictionary",dictionaryObject:"=constructionDictionaryObject",index:"=constructionIndex",editable:"=constructionEditable",fields:"=constructionStep",removeStep:"&constructionRemove",processto:"=constructionProcessto"},compile:function(){return{pre:function(b,f){b.reactions=a.reactions,b.reactionsArray=_.toArray(b.reactions),b.processStep=function(){b.dict=b.dictionaryObject},b.compileStep=function(){function a(a,b,c){return'<div class="span4" construction-field="'+b+'" construction-type="'+a+'" construction-editable="editable" construction-stepindex="index" construction-dictionary="dict" ng-model="'+c+'" ></div>'}b.processStep();var g=angular.element('<div class="constructionStep clearfix"><div class="arrow"></div><div class="reaction" ng-class="{\'errorBackground\' : !reactions[step.reaction]}"><i ng-class="{\'glyphicon glyphicon-resize-full\' : !showReaction, \'glyphicon glyphicon-resize-small\' : showReaction}" style="cursor : pointer;" ng-click="toggleReaction()"></i> {{ reactions[step.reaction].readable || step.reaction }}<i ng-if="editable" class="pull-right glyphicon glyphicon-remove" ng-click="confirmRemoveStep()" style="cursor : pointer; margin-top: 3px"></i> </div><div ng-show="showReaction"><div class="fields" ng-style="{cursor : (editable ? \'move\' : \'inherit\')}"><stepFields class="row clearfix"></stepFields></div><div class="output" ng-class="{\'errorBackground\' : !dictionaryObject[step.output], \'warningBackground\' : !dictionaryObject[step.output].value}" tooltip="{{ index < processto ? (dictionaryObject[step.output].value | stringEnds) : \'<unprocessed>\' }}">Output: <code contenteditable="{{ editable }}" ng-model="step.output" style="display: inline-block; padding: 0 4px;"></code></div></div>');angular.isDefined(b.fields)?angular.forEach(b.fields,function(b){b.model="step."+b.model;var c=a(b.type,b.name,b.model);g.find("stepFields").append(c)}):d.get("views/_dna/construction/steps/"+b.step.reaction+".html",{cache:e}).then(function(a){var d=c(a.data)(b);g.find("stepFields").append(d)},function(){console.log("template not found"),angular.forEach(b.step.input,function(b,c){var d;d=angular.isArray(b)?a("array","","step.input["+c+"]"):"true"==b||"false"==b?a("boolean","","step.input["+c+"]"):a("value","","step.input["+c+"]"),g.find("stepFields").append(d)})}),f.html(c(g)(b))},b.compileStep()},post:function(a){a.$watch("dictionaryObject",function(){a.processStep()}),a.showReaction=!0,a.toggleReaction=function(){a.showReaction=!a.showReaction},a.confirmRemoveStep=function(){var b=h.messageBox("Confirm Remove","Are you sure you want to remove this reaction?",[{label:"Cancel",cssClass:"",result:!1},{label:"Remove",cssClass:"btn-danger",result:!0}]);b.open().then(function(b){b&&a.removeStep()})}}}}}}]),angular.module("clotho.dna").directive("constructionFull",["Construction",function(a){return{restrict:"A",require:"ngModel",replace:!0,scope:{file:"=ngModel",editable:"=constructionEditable",onChange:"&?constructionOnchange",hideopts:"=?constructionHideopts"},templateUrl:"views/_dna/construction/full.html",compile:function(){return{pre:function(){},post:function(b){var c={dictionary:!1,steps:!1,addStep:!1,uptostep:99,processto:99};b.hideopts=angular.extend(c,b.hideopts),b.$watch("file",function(c){c&&a.process(c,!0).then(function(a){console.log("finalResult",a),b.file.dictionary=a,b.onChange({file:b.file})})},!0),b.reactions=a.reactions,b.reactionsArray=_.toArray(b.reactions),b.newStep={reaction:"",input:[],output:""},b.stepValid=function(){return!!a.reactions[b.newStep.reaction]},b.addStep=function(){b.stepValid(b.newStep)&&(b.newStep.input=a.reactions[b.newStep.reaction].input,b.file.steps.push(b.newStep),b.newStep={},b.hideopts.uptostep=b.file.steps.length)}}}}}}]);