angular.module('clotho.editor')
.controller('Editor_FunctionCtrl', function($scope, Clotho, $filter, $q, codemirrorLoader, ClothoSchemas, $timeout) {

	/* data types */

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
			category : 'Primitive',
			reference : false
		});
	});

	ClothoSchemas.retrievedSchemas.then(function (schemas) {
		angular.forEach(schemas, function(schema){
			$scope.paramTypes.push(angular.extend(schema, {category:'Schema'}));
		});
	});

	// todo - this will likely not be practical after release
	$scope.clothoFunctions = [];
	Clotho.query({schema : ClothoSchemas.sharableTypes.Function.schema}, {mute : true}).then(function(result) {
		$scope.clothoFunctions = result;
	});

	$scope.querySchemaWrapper = function(schemaType, value) {
		return Clotho.autocomplete(value).then(function (results) {
			return _.filter(results, function (result) {
				return result.schema == schemaType;
			});
		});
	};

	/* args + deps */

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
		$scope.newDependencyModel = '';
	};

	/* tests */

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

		//fixme - need to resolve arguments which are Clotho objects (they should be IDs)
		//future - will need to do the same for parameterized arguments

		var resolvedArgs = [];
		_.forEach(data.args, function (arg, index) {
			if ($scope.simpleTypes[$scope.sharable.args[index].type]) {
				resolvedArgs.push(arg);
			} else {
				resolvedArgs.push(Clotho.get(arg, {mute : true}).then(function (r) { return r}))
			}
		});

		$q.all(resolvedArgs).then(function (resolved) {
			Clotho.run(data.id, resolved)
			.then(function onFunctionTestSuccess (result){
				$scope.testResults[index] = angular.equals(result, $scope.sharable.tests[index].output.value);
			}, function onFunctionTestError () {
				Clotho.alert('test errored! See activity log.');
			});
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

	//overwrite save to reset tests and set code pristine
	$scope.save = function() {
		$scope.resetTests();
		$scope.Code.$setPristine();
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