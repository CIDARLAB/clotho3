angular.module("clotho.construction",["clotho.foundation","clotho.dna"]),angular.module("clotho.construction").value("ConstructionReactions",{pcr:{reactionId:"clotho.functions.dna.pcr",readable:"PCR",template_step:"views/_construction/step/pcr.html",template_wetlabL:"PCR {{step.input[1]}}/{{step.input[2]}} on {{step.input[0]}}",template_wetlabR:"({{step.output}})"},ligate:{reactionId:"clotho.functions.dna.ligate",readable:"Ligate",template_step:"views/_construction/step/ligate.html",template_wetlabL:"Ligate {{step.input[0] | joinArray:'/'}}",template_wetlabR:"({{step.output}})"},digest:{reactionId:"clotho.functions.dna.digest",readable:"Digest",template_step:"views/_construction/step/digest.html",template_wetlabL:"Digest {{step.input[0]}}",template_wetlabR:"({{step.input[1] | joinArray:'/'}}, {{step.output}})"},gelpurify:{reactionId:"clotho.functions.dna.gelPurify",readable:"Gel Purify",template_step:"views/_construction/step/gelpurify.html",template_wetlabL:"GelPurify {{step.input[0]}}",template_wetlabR:"({{step.input[1] || 'L'}}, {{step.output}})"}}),angular.module("clotho.construction").service("ConstructionSimulator",["ConstructionReactions","$q","Debug","Clotho",function(a,b,c,d){function e(a){return _.isString(a)?d.query({name:a},{mute:!0}).then(function(b){return 1!=b.length&&j.log("ambiguous query! "+b.length+" results for "+a),b[0]}):b.when(a)}function f(a,c){if(!_.isString(c))return b.when(c);if(!_.isEmpty(a.dictionary)){var d=a.dictionary[c];if(!_.isUndefined(d))return b.when(d)}return e(c)}function g(a,b){var c=_.flatten(b.input),d=_.map(c,function(b){return f(a,b)});return _.zipObject(c,d)}function h(a,b){function c(b){return _.map(b,function(b){return _.isString(b)?a.dictionary[b]:_.isArray(b)?c(b):b})}return c(b.input)}function i(b,c){var e=a[c.reaction],f=h(b,c);return d.run(e.reactionId,f,{mute:!0}).then(function(a){var b={};return b[c.output]=a,b})}var j=new c("ConstructionSimulator","#ee8888"),k="final",l=function(a,c){var d=_.assign({dictionary:{}},c?_.clone(a,!0):a);if(_.isEmpty(d.steps))return b.when(d);var f=_.keys(d.dictionary),g={};return j.groupCollapsed("Processing Input Arguments"),_.forEach(d.steps,function(a){f.push(a.output),_.forEach(_.flatten(a.input),function(a){_.isString(a)&&_.indexOf(f,a)<0&&(g[a]=e(a))})}),b.all(g).then(function(a){return j.groupEnd(),_.assign(d.dictionary,a),d})},m=function(a,c){var d=_.assign({dictionary:{}},c?_.clone(a,!0):a);if(_.isEmpty(d.steps))return b.when(d);var e=b.when();return j.groupCollapsed("Simulating File"),_.forEach(d.steps,function(a,c){e=e.then(function(){return j.log("step "+c,d.dictionary,a.input,a.output),b.all(g(d,a))}).then(function(b){return j.log("step "+c+" args resolved",b),_.forEach(b,function(a,b){(_.isArray(a)||_.isObject(a)&&_.isUndefined(d.dictionary[b]))&&(d.dictionary[b]=a)}),i(d,a)}).then(function(a){return j.log("step "+c+" processed",a),b.when(_.assign(d.dictionary,a))})}),e.then(function(){j.groupEnd();var a=_.last(d.steps).output;return a!=k&&(d.dictionary[k]=d.dictionary[a]),_.isArray(d.dictionary[k])?d.dictionary[k][0].sequence=d.dictionary[k][0].sequence.replace(/[\^\|_]/gi,""):_.isObject(d.dictionary[k])&&!_.isUndefined(d.dictionary[k].sequence)&&(d.dictionary[k].sequence=d.dictionary[k].sequence.replace(/[\^\|_]/gi,"")),d})};return{interpolateInputs:l,process:m}}]),angular.module("clotho.construction").directive("constructionFile",["ConstructionSimulator",function(a){return{restrict:"A",templateUrl:"views/_construction/constructionFile.html",scope:{inputFile:"=constructionFile",noProcess:"="},link:function(b,c,d){var e=angular.isUndefined(d.noPreprocess),f=angular.isDefined(d.autoProcess);b.$watch("inputFile",function(a){b.file=a,f?b.process():e&&b.preprocess()}),b.preprocess=function(){a.interpolateInputs(b.file)},b.process=function(){b.processing=!0,a.process(b.file).then(function(){b.processed=!0},function(){b.processError=!0}).finally(function(){b.processing=!1})}}}}]),angular.module("clotho.construction").directive("constructionStep",["ConstructionReactions","$parse",function(a){return{restrict:"A",templateUrl:"views/_construction/constructionStep.html",link:function(b){b.$watch("step",function(c){b.reaction=a[c.reaction]}),b.$watch("file.dictionary[step.output]",function(a){b.stepComputed=angular.isDefined(a)})}}}]).directive("constructionStepInner",["$http","$compile","ConstructionReactions",function(a,b){return{restrict:"A",link:function(c,d){c.$watch("reaction.template_step",function(e){angular.isDefined(e)&&a.get(e,{cache:!0}).then(function(a){d.html(b(a.data)(c))})})}}}]),angular.module("clotho.construction").directive("constructionParserWetlab",function(){return{restrict:"A",scope:{file:"=constructionParserWetlab"},templateUrl:"views/_construction/constructionParserWetlab.html",link:function(){}}}).directive("constructionParserWetlabStep",["ConstructionReactions","$compile",function(a,b){return{restrict:"A",link:function(c,d,e){e.$observe("stepType",function(e){d.html(b("<td>"+a[e].template_wetlabL+"</td><td>"+a[e].template_wetlabR+"</td>")(c))})}}}]);