angular.module('clotho.editor').controller('Editor_SchemaCtrl', function($scope, Clotho) {

	$scope.schemas = [];
	Clotho.query({"schema": "Schema"}).then(function(data) {
		$scope.schemas = data;
	});

	$scope.clothoFunctions = [];
	Clotho.query({"schema": "Function"}).then(function(data) {
		$scope.clothoFunctions = data;
	});

	$scope.accessTypes = [
		{name:'Public', value:'PUBLIC'},
		{name:'Private', value:'PRIVATE'},
		{name:'Read Only', value:'READONLY'}
	];

	$scope.constraintTypes = [
		{name:'RegExp', value:'regex'},
		{name: 'Not Null', value: 'notnull'}
	];

	$scope.primitiveToJava = {
		"string" : "java.lang.String",
		"number" : "java.lang.Long",
		"boolean" : "java.lang.Boolean",
		"object" : "java.util.HashMap",
		"array" : "java.util.List"
	};

	$scope.findSpacesRegExp = /\s/ig;

	$scope.parseField = function(field) {
		//only passed field.value so model maps onto options properly in html
		if ($scope.simpleTypes[field.type]) {
			field.javaType = $scope.primitiveToJava[field.type];
			field.reference = false;
		} else {

			field.reference = true;
		}
	};



	$scope.$watch('sharable.superClass', function() {
		$scope.getSuperClass();
	});

	$scope.getSuperClass = function () {
		if ($scope.sharable.superClass) {
			Clotho.get($scope.sharable.superClass)
				.then(function(result) {
					$scope.superClassObj = result
				})
		}
	};




	$scope.newMethod = function() {
		return ""
	};

	$scope.addMethod = function(method) {
		if (angular.isEmpty($scope.sharable.methods)) {$scope.sharable.methods = [];}
		$scope.sharable.methods.push(method);
	};

	$scope.addNewMethod = function() {
		if (angular.isEmpty($scope.newMethodObj)) return;

		$scope.addMethod($scope.newMethodObj);
		$scope.newMethodObj = $scope.newMethod();
	};

	$scope.determineMethodName = function(id) {
		return _.find($scope.clothoFunctions, {id : id}).name;
	};

	$scope.newField = function() {
		return {
			name: "",
			type: "",
			description: "",
			example: "",
			constraints: null,
			access: "PUBLIC"
		}
	};

	$scope.addField = function() {
		if (angular.isEmpty($scope.sharable.fields)) {$scope.sharable.fields = [];}
		$scope.sharable.fields.push($scope.newField());
	};

	//constraints are processed in saveSchema (link function)
	$scope.newConstraint = function() {
		return {
			type: "",
			value: ""
		};
	};

	$scope.addConstraint = function(index) {
		if (angular.isEmpty($scope.sharable.fields[index].constraints))
			$scope.sharable.fields[index].constraints = [];
		$scope.sharable.fields[index].constraints.push($scope.newConstraint())
	};

	//overwrite the save function for schemas specifically
	$scope.save = function () {
		//process constraints to object
		angular.forEach($scope.sharable.fields, function(field) {
			console.log(field);
			if (field.constraints) {
				var constraintsArray = field.constraints;
				field.constraints = {};
				angular.forEach(constraintsArray, function(constraint) {
					field.constraints[constraint.type] = constraint.value;
				});
			} else {
				field.constraints = null;
			}
		});

		console.log($scope.sharable);

		//inherited
		//verify
		$scope.$parent.save();
	};

});