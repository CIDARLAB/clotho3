angular.module("clotho.editor",["clotho.foundation","clotho.interface"]),angular.module("clotho.editor").directive("clothoEditor",["Clotho","$compile","$parse","$http","$templateCache","$filter","$q","Debug","ClothoSchemas",function(a,b,c,d,e,f,g,h,i){var j=new h("Editor","#dd44dd"),k="views/_editor/_baseFields.html",l="views/_editor/_formActions.html";return{restrict:"A",require:"^form",scope:{sharable:"=?",sharableId:"=?",editMode:"=?",formHorizontal:"=?"},link:function(c,f,g,h){function m(){f.html('<p class="text-center">Please select an Object</p>')}function n(a){return angular.isUndefined(a)||!angular.isObject(a)||angular.isEmpty(a)?(j.warn("editor must be created from object, was given: ",a),void m()):void i.determineSharableType(a).then(function(a){return a},function(){return i.dirtyDetermineType(a)}).then(function(a){var g;console.warn("sharable type is",a),g="Instance"==a?i.sharableTypes.Instance.editor_template_url:i.sharableTypes[a].editor_template_url,d.get(g,{cache:e}).success(function(d){c.showJsonEditor=!1,c.panelClass=i.sharableTypes[a].class||"default";var e=b(d)(c);f.html(e)}).error(function(){j.error("Could not retrieve template: "+g),f.html('<p class="text-center">Please select an Object</p>')})})}function o(b){p(),p=a.watch(b,function(a){angular.equals(c.sharable,a)||(c.sharable=a)},c)}if(angular.isUndefined(g.sharable)&&angular.isUndefined(g.sharableId))return void j.warn("editor not compiling - did not pass either sharable or sharable-id");g.$set("novalidate","novalidate"),angular.isString(g.name)||g.$set("name","sharableEditor"),c.$watch("formHorizontal",function(a){f.toggleClass("form-horizontal",a)}),c._objBaseFields=k,c._formActions=l,i.retrievedSchemas.then(function(a){c.clothoSchemas=a}),c.isValidJson=function(a){var b=!0;try{angular.fromJson(a)}catch(c){b=!1}return b},c.edit=function(){c.editMode=!0},c.toggleJsonEdit=function(){c.showJsonEditor=!c.showJsonEditor},c.reset=function(){a.get(c.sharable.id,{mute:!0}).then(function(a){c.sharable=a,h.$setPristine()})},c.save=function(){c.isValidJson(c.sharable)?(a.set(c.sharable).then(function(a){c.sharable.id=a}),c.editMode=!1):a.alert("Your JSON is invalid... please fix it")},c.discard=function(){c.reset(),c.editMode=!1},c.destroy=function(){a.destroy(c.sharable.id).then(function(){c.editMode=!1,c.sharable=null,m()})},a.listen("collector_reset",function(){m()},c);var p=angular.noop;c.$watch("sharableId",function(b,d){j.log("new id input",b,d),b&&a.get(b,{mute:!0}).then(function(a){c.sharable=a})}),c.$watch("sharable.schema",function(a,b){j.log("sharable schema has changed",a,b),n(c.sharable)}),c.$watch("sharable.id",function(a,b){j.log("sharable id has changed",a,b),o(a),!angular.isEmpty(c.sharable)&&angular.isUndefined(c.sharable.schema)&&n(c.sharable)})}}}]),angular.module("clotho.editor").service("codemirrorLoader",["Clotho","$q",function(a,b){function c(a){var c=f[angular.lowercase(a)];return c?$clotho.extensions.mixin(c):(console.log("language not supported in codemirror: "+a),b.when())}function d(){return b.all([$clotho.extensions.mixin("bower_components/codemirror/lib/codemirror.js"),$clotho.extensions.css("bower_components/codemirror/lib/codemirror.css")])}var e="bower_components/codemirror/mode/",f={javascript:e+"javascript/javascript.js",css:e+"css/css.js",html:e+"htmlmixed/htmlmixed.js",python:e+"python/python.js",groovy:e+"groovy/groovy.js",java:e+"clike/clike.js",markdown:e+"markdown/markdown.js"};return{loadMain:d,loadLanguage:c}}]),angular.module("clotho.editor").directive("formFieldEnumeration",["Clotho","ClothoSchemas","$q","$compile",function(a,b,c,d){function e(a){var c=angular.element("<div ng-form>");return angular.forEach(a,function(a){var d=angular.element("<form-field>");d.attr({name:a.name,horizontal:"formHorizontal"});var e,f=b.formTypeMap[a.type]||!1;if(f&&f.input)e=angular.element("<input>"),e.attr(f.input),e.attr({"ng-model":"sharable."+a.name});else{e=angular.element('<div class="input-group" ng-init="showTypeahead = false">');var g=angular.element("<textarea>");g.attr({"json-edit":"sharable."+a.name,rows:1,placeholder:"Edit JSON directly, use quotes for strings","ng-if":"!showTypeahead","class":"form-control"});var h=angular.element("<input>");h.attr({"ng-if":"showTypeahead",placeholder:"Select object with autocompletion","class":"form-control"});var i=angular.element('<span class="input-group-btn">');i.append(angular.element('<button class="btn btn-default" type="button" ng-click="showTypeahead = !showTypeahead"><span class="glyphicon glyphicon-refresh"></span></button>')),d.attr("no-styling",!0),e.append(g),e.append(h),e.append(i)}e.attr({placeholder:a.description}),a.required&&e.attr({"ng-required":!0}),d.append(e),c.append(d)}),c}return{restrict:"A",scope:{fields:"=?",sharable:"=?",stripBasicFields:"@?"},controller:["$scope","$element","$attrs",function(){}],link:function(c,f,g){function h(){f.empty(),f.prepend(k),f.prepend(j),f.prepend(l)}var i=angular.isDefined(g.hideWarnings);if(angular.isUndefined(g.fields)&&angular.isUndefined(g.sharable))return void(i||f.html("no schema information passed"));c.formHorizontal=angular.isDefined(g.formHorizontal)||c.$parent.formHorizontal;var j=angular.element("<div>"),k=angular.element("<div>"),l=angular.element("<div>");c.$watch("fields",function(a){a&&l.replaceWith(d(e(a))(c))},!0),c.$watch(function(){return _.keys(c.sharable)},function(f){if(c.sharable&&c.sharable.schema)a.get(c.sharable.schema).then(function(a){b.getSuperclassFields(a).then(function(a){var f=_.map(_.keys(c.sharable),function(a){return{name:a}});c.stripBasicFields&&(_.remove(a,function(a){return angular.isDefined(b.sharableBasicFields[a.name])}),_.remove(f,function(a){return angular.isDefined(b.sharableBasicFields[a.name])})),j=d(e(a))(c),k=d(e(f))(c),h()})});else{j=angular.element(i?"":'<div class="alert alert-warning">Sharable has no schema...</div>');var g=f;_.remove(g,function(a){return angular.isDefined(b.sharableBasicFields[a])});var l=_.map(g,function(a){return{name:a}});k=d(e(l))(c),h()}},!0)}}}]),angular.module("clotho.editor").directive("formFieldInput",["ClothoSchemas","$compile","$parse",function(a,b,c){return{restrict:"A",link:function(d,e,f){function g(c){var g,h=a.formTypeMap[c]||!1;h?h.input?(g=angular.element("<input>"),g.attr(h.input),g.attr({"ng-model":f.formFieldModel})):(g=angular.element("<textarea>"),g.attr({"json-edit":f.formFieldModel,rows:3})):(g=angular.element("<input>"),g.attr({type:"text","ng-model":f.formFieldModel})),angular.forEach(f,function(a,b){"class"==b?g.addClass(a):"$"==b.charAt(0)||g.attr(b)||g.attr(f.$attr[b],a)}),e.attr("type",null),e.attr("ng-model",null),e.attr("json-edit",null),e.replaceWith(b(g)(d))}var h;h=d.$watch(function(){return c(f.formFieldType)(d)},function(a,b){a!=b&&(g(a),h())})}}}]),angular.module("clotho.editor").directive("jsonEdit",function(){return{restrict:"A",require:"ngModel",template:'<textarea ng-model="jsonEditing"></textarea>',replace:!0,scope:{model:"=jsonEdit"},link:function(a,b,c,d){function e(a){try{return angular.fromJson(a)}catch(b){return j(),a}}function f(a){return angular.toJson(a,!0)}function g(b){a.jsonEditing=f(b)}function h(b){a.model=e(b)}function i(){d.$setValidity("json",!0)}function j(){d.$setValidity("json",!1)}function k(a){var b=angular.isDefined(a);try{angular.fromJson(a)}catch(c){b=!1}return b}g(a.model),a.$watch("jsonEditing",function(a){k(a)?(i(),h(a)):j()},!0),a.$watch("model",function(a){g(a)},!0)}}}).directive("jsonEditor",function(){return{restrict:"A",require:"ngModel",link:function(a,b,c,d){function e(a){try{var b=angular.fromJson(a);return d.$setValidity("json",!0),b}catch(c){return d.$setValidity("json",!1),a}}function f(a){return angular.toJson(a,!0)}d.$parsers.push(e),d.$formatters.push(f)}}}),angular.module("clotho.editor").controller("Editor_SharableCtrl",["$scope",function(a){a.addNewField=function(){a.newFieldKey&&a.newFieldVal&&(a.sharable[a.newFieldKey]=a.newFieldVal,a.newFieldKey=null,a.newFieldVal=null)},a.checkNewFieldName=function(){angular.isDefined(a.sharable[a.newFieldKey])&&a.newField.$setInvalid()}}]).directive("sharableCheckFieldExists",["$parse",function(a){return{restrict:"A",require:"ngModel",link:function(b,c,d,e){e.$parsers.push(function(c){if(c){var f=a(d.sharableCheckFieldExists)(b),g=angular.isDefined(f[c]);return e.$setValidity("duplicateField",!g),c}})}}}]),angular.module("clotho.editor").controller("Editor_FunctionCtrl",["$scope","Clotho","$filter","$q","codemirrorLoader","ClothoSchemas","$timeout",function(a,b,c,d,e,f,g){a.langTypes=[{name:"JavaScript",value:"JAVASCRIPT"},{name:"Python",value:"PYTHON"}],a.outputTypes=[{name:"Value",value:"VALUE"},{name:"Reference",value:"REFERENCE"}],a.simpleTypes={object:!0,string:!0,number:!0,"boolean":!0,array:!0},a.paramTypes=[],angular.forEach(f.primitiveToJava,function(b,c){a.paramTypes.push({id:c,name:c,type:c,category:"Primitive",reference:!1})}),f.retrievedSchemas.then(function(b){angular.forEach(b,function(b){a.paramTypes.push(angular.extend(b,{category:"Schema"}))})}),a.clothoFunctions=[],b.query({schema:f.sharableTypes.Function.schema},{mute:!0}).then(function(b){a.clothoFunctions=b}),a.querySchemaWrapper=function(a,c){return b.autocomplete(c).then(function(b){return _.filter(b,function(b){return b.schema==a})})},a.addArg=function(){angular.isEmpty(a.sharable.args)&&(a.sharable.args=[]),a.sharable.args.push({type:"",name:""})},a.addDep=function(){angular.isEmpty(a.sharable.dependencies)&&(a.sharable.dependencies=[]),""!=a.newDependencyModel&&a.sharable.dependencies.push(a.newDependencyModel),a.newDependencyModel=""},a.addTest=function(){angular.isEmpty(a.sharable.tests)&&(a.sharable.tests=[]),a.sharable.tests.push({args:[],output:{value:"",type:""}})},a.testResults={},a.singleTest=function(c){var e={};e.id=a.sharable.id,e.args=angular.isEmpty(a.sharable.tests)?[]:a.sharable.tests[c].args;var f=[];_.forEach(e.args,function(c,d){f.push(a.simpleTypes[a.sharable.args[d].type]?c:b.get(c,{mute:!0}).then(function(a){return a}))}),d.all(f).then(function(d){b.run(e.id,d).then(function(b){a.testResults[c]=angular.equals(b,a.sharable.tests[c].output.value)},function(){b.alert("test errored! See activity log.")})})},a.runAllTests=function(){for(var b=0;b<a.sharable.tests.length;b++)a.singleTest(b)},a.resetTests=function(){a.testResults={}},a.save=function(){a.resetTests(),a.Code.$setPristine(),a.$parent.save()},a.$watch("editMode",function(b){a.codemirrorEditorOptions.readOnly=b?!1:"nocursor",a.resetTests()}),a.codemirrorEditorOptions={lineWrapping:!0,lineNumbers:!0,onLoad:function(b){a.$watch("sharable.language",function(a){if(a){var c=a.toLowerCase();e.loadLanguage(c).then(function(){c="java"==c?"text/x-java":c,g(function(){b.setOption("mode",c)},500)})}}),b.on("beforeChange",function(){}),b.on("change",function(){})}},a.resetTests()}]),angular.module("clotho.editor").controller("Editor_TrailCtrl",function(){}),angular.module("clotho.editor").controller("Editor_SchemaCtrl",["$scope","Clotho","ClothoSchemas",function(a,b,c){a.schemas=[],c.retrievedSchemas.then(function(b){a.schemas=b}),a.clothoFunctionWrap=function(a){return b.query({schema:c.sharableTypes.Function.schema,name:a},{mute:!0}).then(function(a){return a||[]})},a.accessTypes=c.accessTypes,a.constraintTypes=c.constraintTypes,a.primitiveToJava=c.primitiveToJava,a.findSpacesRegExp=/\s/gi,a.paramTypes=[{name:"object",type:"object",category:"Primitive",reference:!1},{name:"array",type:"array",category:"Primitive",reference:!1},{name:"string",type:"string",category:"Primitive",reference:!1},{name:"number",type:"number",category:"Primitive",reference:!1},{name:"boolean",type:"boolean",category:"Primitive",reference:!1}],c.retrievedSchemas.then(function(b){angular.forEach(b,function(b){a.paramTypes.push(angular.extend(b,{category:"Schema"}))})}),a.$watch("sharable.superClass",function(b){b&&a.getSuperClass()}),a.getSuperClass=function(){c.getSuperclassFields(a.sharable).then(function(b){a.superClassFields=b,a.superClassObj=b})},a.newMethod=function(){return""},a.addMethod=function(b){angular.isEmpty(a.sharable.methods)&&(a.sharable.methods=[]),a.sharable.methods.push(b)},a.addNewMethod=function(){angular.isEmpty(a.newMethodObj)||(a.addMethod(a.newMethodObj),a.newMethodObj=a.newMethod())},a.determineMethodName=function(b){return _.find(a.clothoFunctions,{id:b}).name},a.newField=function(){return{name:"",type:"",description:"",example:"",constraints:null,access:"PUBLIC"}},a.addField=function(){angular.isEmpty(a.sharable.fields)&&(a.sharable.fields=[]),a.sharable.fields.push(a.newField())},a.newConstraint=function(){return{type:"",value:""}},a.addConstraint=function(b){angular.isEmpty(a.sharable.fields[b].constraints)&&(a.sharable.fields[b].constraints=[]),a.sharable.fields[b].constraints.push(a.newConstraint())},a.save=function(){angular.forEach(a.sharable.fields,function(a){if(console.log(a),a.constraints){var b=a.constraints;a.constraints={},angular.forEach(b,function(b){a.constraints[b.type]=b.value})}else a.constraints=null}),a.$parent.save()}}]);