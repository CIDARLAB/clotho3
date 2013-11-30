//todo - handle loading of markdown, codemirror, etc.
//todo - handle initial coloring of code based on language

angular.module('clotho.editor').controller('Editor_FunctionCtrl', function($scope, Clotho, $filter) {

	$scope.langTypes = [
		{name:'JavaScript', value:'JAVASCRIPT'},
		{name:'Java', value:'JAVA'},
		{name:'Python', value:'PYTHON'},
		{name:'Groovy', value:'GROOVY'}
	];

	$scope.outputTypes = [
		{name:'Value', value:'VALUE'},
		{name:'Reference', value:'REFERENCE'}
	];

	$scope.simpleTypes = {
		"object" : true,
		"string" : true,
		"number" : true,
		"boolean" : true,
		"array" : true
	};

	$scope.paramTypes = [
		{name:'object', type: "object", category:'Primitive', javaType : "java.util.HashMap", reference: false},
		{name:'array', type : "array", category:'Primitive', javaType : "java.util.Arrays", reference: false},
		{name:'string', type : "string", category:'Primitive', javaType : "java.lang.String", reference: false},
		{name:'number', type : "number", category:'Primitive', javaType : "java.lang.Long", reference: false},
		{name:'boolean', type : "boolean", category:'Primitive', javaType : "java.lang.Boolean", reference: false}
	];
	Clotho.query({schema:"Schema"}).then(function(data){
		angular.forEach(data, function(schema){
			$scope.paramTypes.push(angular.extend(schema, {category:'Schema'}));
		});
	});

	$scope.clothoFunctions = [];
	Clotho.query({schema: "Function"}).then(function(result) {
		$scope.clothoFunctions = result;
	});



	$scope.addArg = function() {
		if (angular.isEmpty($scope.editable.args)) {$scope.editable.args = [];}
		$scope.editable.args.push({"type" : "", "name" : ""});
	};

	$scope.addDep = function() {
		if (angular.isEmpty($scope.editable.dependencies)) {$scope.editable.dependencies = [];}
		$scope.editable.dependencies.push("");
	};

	$scope.addTest = function() {
		if (angular.isEmpty($scope.editable.tests)) {$scope.editable.tests = [];}
		$scope.editable.tests.push({"args" : [], "output" : {"value" : "", "type" : ""}});
	};

	$scope.testResults = {};
	$scope.singleTest = function(index) {

		var data = {};
		data.id = $scope.editable.id;
		if (angular.isEmpty($scope.editable.tests)) {
			data.args = [];
		}
		else {
			data.args = $scope.editable.tests[index].args;
		}

		Clotho.run(data.id, data.args).then(function (result){
			console.log(result, $scope.editable.tests[index].output.value, result == $scope.editable.tests[index].output.value);
			$scope.testResults[index] = (result == $scope.editable.tests[index].output.value);
			/* if (result == angular.fromJson($scope.editable.testResult)) {
			 ClientAPI.say({text:"test success!"});
			 } else {
			 ClientAPI.say({text:"test failed!"});
			 }*/
		});
	};

	$scope.runAllTests = function() {
		for (var i = 0; i < $scope.editable.tests.length; i++) {
			$scope.singleTest(i);
		}
	};

	$scope.resetTests = function() {
		$scope.testResults = {};
	};

	$scope.querySchemaWrapper = function(schemaType) {
		return Clotho.query({schema: schemaType}).then(function (result) {
			return $filter('limitTo')(result, 10);
		})
	};


	// code mirror

	$scope.$watch('editMode', function(newval) {
		$scope.codemirrorEditorOptions.readOnly = (newval) ? false : 'nocursor';
	});

	$scope.codemirrorEditorOptions = {
		lineWrapping : true,
		lineNumbers: true,
		onLoad : function(_cm){
			// HACK to have the codemirror instance in the scope...
			$scope.$watch('editable.language', function(newlang) {
				// HACK to catch java case
				var mode = $scope.editable.language.toLowerCase();
				mode = (mode == 'java') ? 'text/x-java' : mode;

				_cm.setOption("mode", mode);
			});

			// example Events
			_cm.on("beforeChange", function(){});
			_cm.on("change", function(){});
		}
	};

	//init()
	$scope.resetTests();
});