angular.module('clotho.editor').controller('Editor_SchemaCtrl', function($scope, Clotho, ClothoSchemas) {

    $scope.schemas = [];

    ClothoSchemas.retrievedSchemas.then(function (schemas) {
        $scope.schemas = schemas;
    });

    $scope.clothoFunctionWrap = function (query) {
        return Clotho.query({"schema": ClothoSchemas.sharableTypes.Function.schema, name : query}, {mute : true})
        .then(function(data) {
            return data || [];
        });
    };

    $scope.accessTypes = ClothoSchemas.accessTypes;

    $scope.constraintTypes = ClothoSchemas.constraintTypes;

    $scope.primitiveToJava = ClothoSchemas.primitiveToJava;

    $scope.findSpacesRegExp = /\s/ig;

    //TODO - dry. duplicates function editor code
    $scope.paramTypes = [
        {name:'object', type: "object", category:'Primitive', reference: false},
        {name:'array', type : "array", category:'Primitive', reference: false},
        {name:'string', type : "string", category:'Primitive', reference: false},
        {name:'number', type : "number", category:'Primitive', reference: false},
        {name:'boolean', type : "boolean", category:'Primitive', reference: false}
    ];

    ClothoSchemas.retrievedSchemas.then(function (schemas) {
        angular.forEach(schemas, function(schema){
            $scope.paramTypes.push(angular.extend(schema, {category:'Schema'}));
        });
    });


    $scope.$watch('sharable.superClass', function(val) {
        if (!!val) {
            $scope.getSuperClass();
        }
    });

    $scope.getSuperClass = function () {
        ClothoSchemas.getSuperclassFields($scope.sharable)
        .then(function(result) {
            $scope.superClassFields = result;
            $scope.superClassObj = result;
        });
    };

    // note - currently unused as methods section in schema.html commented out
    // $scope.newMethod = function() {
    //     return ""
    // };

    // $scope.addMethod = function(method) {
    //     if (angular.isEmpty($scope.sharable.methods)) {$scope.sharable.methods = [];}
    //     $scope.sharable.methods.push(method);
    // };

    // $scope.addNewMethod = function() {
    //     if (angular.isEmpty($scope.newMethodObj)) return;

    //     $scope.addMethod($scope.newMethodObj);
    //     $scope.newMethodObj = $scope.newMethod();
    // };

    // $scope.determineMethodName = function(id) {
    //  return _.find($scope.clothoFunctions, {id : id}).name;
    // };

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
        $scope.$parent.save();
    };

});
