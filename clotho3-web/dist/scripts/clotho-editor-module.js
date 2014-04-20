angular.module("clotho.editor",["clotho.foundation","clotho.interface"]),angular.module("clotho.editor").directive("clothoEditor",["Clotho","$compile","$parse","$http","$templateCache","$filter","$q","Debug","ClothoSchemas",function(a,b,c,d,e,f,g,h,i){var j=new h("Editor","#dd44dd");return{restrict:"A",require:["ngModel","^form"],scope:{inputSharable:"=ngModel",editMode:"=?"},controller:["$scope","$element","$attrs",function(c,f,h){return h.$set("novalidate","novalidate"),angular.isString(h.name)||h.$set("name","sharableEditor"),c.logScope=function(){j.log(c)},c._objBaseFields="views/_editor/_baseFields.html",c._formActions="views/_editor/_formActions.html",c.getPartialAndCompile=function(a){return c.showJsonEditor=!1,angular.isEmpty(a)?void f.html('<p class="text-center">Please select an Object</p>'):void d.get(i.sharableTypes[a].editor_template_url,{cache:e}).success(function(a){var d=b(a)(c);f.html(d)}).error(function(){j.error("Could not retrieve template for type "+a),f.html('<p class="text-center">Please select an Object</p>')})},c.processInputSharable=function(b){if(c.sharable={},angular.isEmpty(b))return j.log("sharable is undefined / empty"),void c.getPartialAndCompile();j.debug("processing sharable: ",b),angular.isObject(b)&&!angular.isEmpty(b)?(j.log("sharable object passed"),c.id=""|b.id,c.sharable=b):angular.isString(b)?(j.log("passed a string, assuming a valid ID: "+b),c.id=b):(j.warn("sharable must be an object or a string, you passed: ",b),c.id="",c.sharable={});var d=angular.isEmpty(c.sharable)&&""!=c.id?a.get(c.id):g.when(c.sharable);d.then(function(a){c.sharable=a,c.getPartialAndCompile(i.determineType(a),a)})},c.processInputSharable(c.inputSharable),{removeField:function(a){j.warn("(DEPRECATED) - removing fields from sharables is not supported"),j.log("removing key "+a),delete c.sharable[a]}}}],compile:function(){return{pre:function(){},post:function(b,c,d,e){b.isValidJson=function(a){var b=!0;try{angular.fromJson(a)}catch(c){b=!1}return b},b.formCtrl=e[1],b.edit=function(){b.editMode=!0},b.toggleJsonEdit=function(){b.showJsonEditor=!b.showJsonEditor},b.reset=function(){b.formCtrl.$setPristine(),a.get(b.id).then(function(a){b.sharable=a})},b.save=function(){b.isValidJson(b.sharable)?(a.set(b.sharable),b.editMode=!1):a.alert("Your JSON is invalid... please fix it")},b.discard=function(){b.reset(),b.editMode=!1},b.destroy=function(){a.destroy(b.sharable.id).then(function(){b.editMode=!1,b.sharable=null,b.id=null,b.processInputSharable({})})},a.listen("collector_reset",function(){b.processInputSharable()},b),a.watch(b.id,b.sharable,b),b.$watch("inputSharable",function(a,c){(!c||a&&c&&a.id!=c.id)&&(b.editMode=!1),b.processInputSharable(a)}),b.$watch("editMode",function(a){j.log("edit mode: ",a)})}}}}}]),angular.module("clotho.editor").service("codemirrorLoader",["Clotho","$q",function(a,b){function c(a){var c=f[angular.lowercase(a)];return c?$clotho.extensions.mixin(c):(console.log("language not supported in codemirror: "+a),b.when())}function d(){return b.all([$clotho.extensions.mixin("bower_components/codemirror/lib/codemirror.js"),$clotho.extensions.css("bower_components/codemirror/lib/codemirror.css")])}var e="bower_components/codemirror/mode/",f={javascript:e+"javascript/javascript.js",css:e+"css/css.js",html:e+"htmlmixed/htmlmixed.js",python:e+"python/python.js",groovy:e+"groovy/groovy.js",java:e+"clike/clike.js",markdown:e+"markdown/markdown.js"};return{loadMain:d,loadLanguage:c}}]),angular.module("clotho.editor").directive("formField",function(){var a='<div class="form-group" ng-form ng-transclude></div>';return{restrict:"E",template:a,replace:!0,transclude:!0,require:["^form","^?clothoEditor"],controller:["$scope","$element","$attrs",function(){}],link:function(a,b,c,d){var e=d[0],f=(d[1],c.name),g=c.help,h=(c.removableField,c.hideLabel),i=b.children();if(1!==i.length)throw"You must include the form input element as the single child of this directive, instead got "+b.html()+"generating:";var j=i.attr("id"),k=(i.prop("tagName"),i.attr("type")),l=/checkbox|radio/gi,m=/file|checkbox|radio/gi;j||(j="input"+Math.floor(1e7*Math.random()).toString(),i.attr("id",j)),m.test(k)||i.addClass("form-control");var n=angular.element('<label class="control-label" for="'+j+'">'+("radio"===k?i.val():f)+"</label>");if(l.test(k)){n.prepend(i);var o=angular.element('<div class="'+k+'"></div>');b.append(o),h=!1}f&&!h&&b.prepend(n),g&&b.append('<p class="help-block">'+g+"</p>"),a.$watch(function(){return e.$invalid},function(a){b.toggleClass("has-error",a)})}}}),angular.module("clotho.interface").directive("functionCodeDrop",["$upload","$window","$timeout",function(a,b,c){return{restrict:"A",templateUrl:"views/_interface/codeDrop.html",scope:{updateOnRead:"="},controller:["$scope","$element","$attrs",function(a){a.uploadRightAway=!1,a.selectedFiles=[],a.inputFiles=[],a.onFileSelect=function(d){function e(b){b.readAsText(d[f]),b.onload=function(b){c(function(){a.updateOnRead=b.target.result})}}console.log("files changed"),a.selectedFiles=d;for(var f=0;f<d.length;f++){var g=d[f];if(b.FileReader&&g.type.indexOf("text")>-1){var h=new FileReader;e(h,f)}}}}],link:function(){}}}]),angular.module("clotho.editor").directive("formFieldEnumeration",["Clotho","ClothoUtils","$q","$compile",function(a,b,c,d){function e(a){var b="",c={string:"text","boolean":"checkbox",select:"select",radio:"radio",object:"object"};return angular.forEach(a,function(a){var d=c[a.type]||a.type||"text";"?"==d&&"text"==a.type;var e,f=a.required?"required":"",g='<form-field name="'+a.name+'" removable-field="'+a.name+'">',h="</form-field>";switch(d){case"object":e='<textarea json-edit="sharable.'+a.name+'" rows="3" '+f+"></textarea>";break;case"radio":angular.forEach(a.options,function(b){e+='<input type="radio" value="'+b+'" ng-model="sharable.'+a.name+">"+b});break;case"select":var i="";angular.forEach(a.options,function(a){i=i+'<option value="'+a+'">'+a+"</option>"}),e="<select "+f+' ng-model="sharable.'+a.name+'">'+i+"</select>";break;default:e='<input type="'+d+'" '+f+'ng-model="sharable.'+a.name+'" >'}b+=g+e+h}),b}return{restrict:"A",scope:{fields:"=?",schema:"=?",sharable:"=?"},controller:["$scope","$element","$attrs",function(a,b){a.generateFields=function(c){var f=c.fields||a.fields||[],g=c.schema||a.schema||{fields:[]},h=c.sharable||a.sharable||{schema:{fields:[]}},i=f.concat(g.fields);_.forEach(h,function(a,b){i.push({name:b})});var j=_.uniq(i,function(a){return a.name});b.html(e(j)),d(b.contents())(a)}}],link:function(c,d){function e(a){return b.downloadSchemaDependencies(a).then(function(a){c.generateFields({schema:a})})}return angular.isUndefined(c.fields)&&angular.isUndefined(c.schema)&&angular.isUndefined(c.sharable)?void d.html("no information passed"):(c.$watch("fields",function(a){a&&(console.log("updating fields"),c.generateFields({fields:a}))},!0),c.$watch("schema",function(a){a&&a.fields&&a.fields.length&&(console.log("updating schema"),e(a))},!0),void c.$watch(function(){return _.keys(c.sharable)},function(){c.sharable&&(console.log("updating sharable"),a.get(c.sharable.schema).then(function(a){e(a)}))},!0))}}}]),angular.module("clotho.editor").directive("jsonEdit",function(){return{restrict:"A",require:"ngModel",template:'<textarea ng-model="jsonEditing"></textarea>',replace:!0,scope:{model:"=jsonEdit"},link:function(a,b,c,d){function e(b){a.jsonEditing=angular.copy(j(b))}function f(b){a.model=i(b)}function g(){d.$setValidity("json",!0)}function h(){d.$setValidity("json",!1)}function i(a){try{return angular.fromJson(a)}catch(b){return h(),a}}function j(a){return angular.toJson(a,!0)}function k(a){var b=!0;try{angular.fromJson(a)}catch(c){b=!1}return b}e(a.model),a.$watch("jsonEditing",function(a,b){a!=b&&(k(a)?(g(),f(a)):h())},!0),a.$watch("model",function(a,b){a!=b&&e(a)},!0)}}}),angular.module("clotho.editor").controller("Editor_SharableCtrl",["$scope",function(a){a.addNewField=function(){a.newFieldKey&&a.newFieldVal&&(a.sharable[a.newFieldKey]=a.newFieldVal,a.newFieldKey=null,a.newFieldVal=null)}}]),angular.module("clotho.editor").controller("Editor_FunctionCtrl",["$scope","Clotho","$filter","codemirrorLoader","ClothoSchemas",function(a,b,c,d,e){a.langTypes=[{name:"JavaScript",value:"JAVASCRIPT"},{name:"Java",value:"JAVA"},{name:"Python",value:"PYTHON"},{name:"Groovy",value:"GROOVY"}],a.outputTypes=[{name:"Value",value:"VALUE"},{name:"Reference",value:"REFERENCE"}],a.simpleTypes={object:!0,string:!0,number:!0,"boolean":!0,array:!0},a.paramTypes=[{name:"object",type:"object",category:"Primitive",javaType:"java.util.HashMap",reference:!1},{name:"array",type:"array",category:"Primitive",javaType:"java.util.Arrays",reference:!1},{name:"string",type:"string",category:"Primitive",javaType:"java.lang.String",reference:!1},{name:"number",type:"number",category:"Primitive",javaType:"java.lang.Long",reference:!1},{name:"boolean",type:"boolean",category:"Primitive",javaType:"java.lang.Boolean",reference:!1}],e.retrievedSchemas.then(function(b){angular.forEach(b,function(b){a.paramTypes.push(angular.extend(b,{category:"Schema"}))})}),a.clothoFunctions=[],b.query({schema:"Function"}).then(function(b){a.clothoFunctions=b}),a.addArg=function(){angular.isEmpty(a.sharable.args)&&(a.sharable.args=[]),a.sharable.args.push({type:"",name:""})},a.addDep=function(){angular.isEmpty(a.sharable.dependencies)&&(a.sharable.dependencies=[]),a.sharable.dependencies.push("")},a.addTest=function(){angular.isEmpty(a.sharable.tests)&&(a.sharable.tests=[]),a.sharable.tests.push({args:[],output:{value:"",type:""}})},a.testResults={},a.singleTest=function(c){var d={};d.id=a.sharable.id,d.args=angular.isEmpty(a.sharable.tests)?[]:a.sharable.tests[c].args,b.run(d.id,d.args).then(function(b){console.log(b,a.sharable.tests[c].output.value,b==a.sharable.tests[c].output.value),a.testResults[c]=b==a.sharable.tests[c].output.value})},a.runAllTests=function(){for(var b=0;b<a.sharable.tests.length;b++)a.singleTest(b)},a.resetTests=function(){a.testResults={}},a.querySchemaWrapper=function(a){return b.query({schema:a}).then(function(a){return c("limitTo")(a,10)})},a.testPopoverText=function(b){var d=a.sharable.tests[b].args;return d?c("json")(d):"Search for an object of type {{ param.type }}"},a.save=function(){a.resetTests(),a.$parent.save()},a.$watch("editMode",function(b){a.codemirrorEditorOptions.readOnly=b?!1:"nocursor"}),a.codemirrorEditorOptions={lineWrapping:!0,lineNumbers:!0,onLoad:function(b){a.$watch("sharable.language",function(a,c){if(a&&a!=c){var e=a.toLowerCase();d.loadLanguage(e).then(function(){e="java"==e?"text/x-java":e,b.setOption("mode",e)})}}),b.on("beforeChange",function(){}),b.on("change",function(){a.resetTests()})}},a.resetTests()}]),angular.module("clotho.editor").controller("Editor_TrailCtrl",function(){}),angular.module("clotho.editor").controller("Editor_SchemaCtrl",["$scope","Clotho","ClothoSchemas",function(a,b,c){a.schemas=[],c.retrievedSchemas.then(function(b){a.schemas=b}),a.clothoFunctions=[],b.query({schema:"Function"}).then(function(b){a.clothoFunctions=b}),a.accessTypes=c.accessTypes,a.constraintTypes=c.constraintTypes,a.primitiveToJava=c.primitiveToJava,a.findSpacesRegExp=/\s/gi,a.parseField=function(b){a.simpleTypes[b.type]?(b.javaType=a.primitiveToJava[b.type],b.reference=!1):b.reference=!0},a.$watch("sharable.superClass",function(){a.getSuperClass()}),a.getSuperClass=function(){a.sharable.superClass&&b.get(a.sharable.superClass).then(function(b){a.superClassObj=b})},a.newMethod=function(){return""},a.addMethod=function(b){angular.isEmpty(a.sharable.methods)&&(a.sharable.methods=[]),a.sharable.methods.push(b)},a.addNewMethod=function(){angular.isEmpty(a.newMethodObj)||(a.addMethod(a.newMethodObj),a.newMethodObj=a.newMethod())},a.determineMethodName=function(b){return _.find(a.clothoFunctions,{id:b}).name},a.newField=function(){return{name:"",type:"",description:"",example:"",constraints:null,access:"PUBLIC"}},a.addField=function(){angular.isEmpty(a.sharable.fields)&&(a.sharable.fields=[]),a.sharable.fields.push(a.newField())},a.newConstraint=function(){return{type:"",value:""}},a.addConstraint=function(b){angular.isEmpty(a.sharable.fields[b].constraints)&&(a.sharable.fields[b].constraints=[]),a.sharable.fields[b].constraints.push(a.newConstraint())},a.save=function(){angular.forEach(a.sharable.fields,function(a){if(console.log(a),a.constraints){var b=a.constraints;a.constraints={},angular.forEach(b,function(b){a.constraints[b.type]=b.value})}else a.constraints=null}),console.log(a.sharable),a.$parent.save()}}]);