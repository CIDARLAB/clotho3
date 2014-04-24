//todo - handle loading of markdown, codemirror, etc.
//todo - handle initial coloring of code based on language

angular.module('clotho.editor').controller('Editor_FunctionCtrl', function($scope, Clotho, $filter, codemirrorLoader, ClothoSchemas, $timeout) {

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

	$scope.paramTypes = [];

	angular.forEach(ClothoSchemas.primitiveToJava, function (val, key) {
		$scope.paramTypes.push({
			id : key,
			name : key,
			type : key,
			javaType : val,
			category : 'Primitive',
			reference : false
		});
	});

	ClothoSchemas.retrievedSchemas.then(function (schemas) {
		angular.forEach(schemas, function(schema){
			$scope.paramTypes.push(angular.extend(schema, {category:'Schema'}));
			console.log($scope.paramTypes);
		});
	});

	// todo - this will likely not be practical after release
	$scope.clothoFunctions = [];
	Clotho.query(ClothoSchemas.sharableTypes.Function.schema).then(function(result) {
		$scope.clothoFunctions = result;
	});

	$scope.addArg = function() {
		if (angular.isEmpty($scope.sharable.args)) {
			$scope.sharable.args = [];
		}
		$scope.sharable.args.push({"type" : "", "name" : ""});
	};

	$scope.addDep = function() {
		if (angular.isEmpty($scope.sharable.dependencies)) {
			$scope.sharable.dependencies = [];
		}
		if ($scope.newDependencyModel != '') {
			$scope.sharable.dependencies.push($scope.newDependencyModel);
		}
	};

	$scope.addTest = function() {
		if (angular.isEmpty($scope.sharable.tests)) {$scope.sharable.tests = [];}
		$scope.sharable.tests.push({"args" : [], "output" : {"value" : "", "type" : ""}});
	};

	$scope.testResults = {};
	$scope.singleTest = function(index) {

		var data = {};
		data.id = $scope.sharable.id;
		if (angular.isEmpty($scope.sharable.tests)) {
			data.args = [];
		}
		else {
			data.args = $scope.sharable.tests[index].args;
		}

		Clotho.run(data.id, data.args).then(function (result){
			console.log(result, $scope.sharable.tests[index].output.value, result == $scope.sharable.tests[index].output.value);
			$scope.testResults[index] = (result == $scope.sharable.tests[index].output.value);
			/* if (result == angular.fromJson($scope.sharable.testResult)) {
			 ClientAPI.say({text:"test success!"});
			 } else {
			 ClientAPI.say({text:"test failed!"});
			 }*/
		});
	};

	$scope.runAllTests = function() {
		for (var i = 0; i < $scope.sharable.tests.length; i++) {
			$scope.singleTest(i);
		}
	};

	$scope.resetTests = function() {
		$scope.testResults = {};
	};

	$scope.querySchemaWrapper = function(schemaType, value) {
		return Clotho.query({schema: schemaType, name : value}, {maxResults : 10})
		.then(function (results) {
			return results;
		});
	};

	$scope.testPopoverText = function (ind) {
		var args = $scope.sharable.tests[ind].args;
		return !!args ? $filter('json')(args) : 'Search for an object of type {{ param.type }}';
	};

	//todo - re-parse non-simple tests - interpolate strings to objects so run correctly

	//overwrite save to reset tests
	$scope.save = function() {
		$scope.resetTests();
		$scope.$parent.save()
	};


	// code mirror

	$scope.$watch('editMode', function(newval) {
		$scope.codemirrorEditorOptions.readOnly = (newval) ? false : 'nocursor';
		$scope.resetTests();
	});

	$scope.codemirrorEditorOptions = {
		lineWrapping : true,
		lineNumbers: true,
		// HACK to have the codemirror instance in the scope...
		onLoad : function(_cm){
			$scope.$watch('sharable.language', function(newlang, oldlang) {
				if (!!newlang) {
					var mode = newlang.toLowerCase();

					codemirrorLoader.loadLanguage(mode).then(function () {
						// HACK to catch java case
						mode = (mode == 'java') ? 'text/x-java' : mode;
						//give it time to register
						$timeout(function () {
							_cm.setOption("mode", mode);
						}, 500);
					});
				}
			});

			// example Events
			_cm.on("beforeChange", function(){});
			_cm.on("change", function(){});
		}
	};

	//init()
	$scope.resetTests();
});