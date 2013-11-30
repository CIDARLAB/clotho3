angular.module('clotho.editor').controller('Editor_SchemaCtrl', function($scope, Clotho) {

	$scope.schemas = [];
	Clotho.query({"schema": "Schema"}).then(function(data) {
		$scope.schemas = data;
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

	$scope.newMethod = function() {
		return ""
	};

	$scope.addMethod = function(method) {
		if (angular.isEmpty($scope.editable.methods)) {$scope.editable.methods = [];}
		$scope.editable.methods.push(method);
	};

	$scope.addNewMethod = function() {
		if (angular.isEmpty($scope.newMethodObj)) return;

		$scope.addMethod($scope.newMethodObj);
		$scope.newMethodObj = $scope.newMethod();
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
		if (angular.isEmpty($scope.editable.fields)) {$scope.editable.fields = [];}
		$scope.editable.fields.push($scope.newField());
	};

	//note - constraints are processed in saveSchema (link function)
	$scope.newConstraint = function() {
		return {
			type: "",
			value: ""
		};
	};

	$scope.addConstraint = function(index) {
		if (angular.isEmpty($scope.editable.fields[index].constraints))
			$scope.editable.fields[index].constraints = [];
		$scope.editable.fields[index].constraints.push($scope.newConstraint())
	};

	$scope.saveSchema = function () {
		//process constraints to object
		angular.forEach($scope.editable.fields, function(field) {
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

		console.log($scope.editable);

		//inherited
		$scope.save();
	};

});