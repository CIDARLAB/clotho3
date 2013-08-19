'use strict';

Application.Schemas.controller('SchemasCtrl', ['$scope', 'Clotho', '$filter', '$location', function($scope, Clotho, $filter, $location) {

//init

    $scope.nav = ["Browse", "Create"];

    $scope.setCurrent = function(mode) {
        $scope.current = mode;
    };


    $scope.newDescription = function(){
        return  {
            name: "NewSchema",
            parentSchema: null,
            description: "",
            fields: []
        };
    };

    $scope.newFieldDescription = function (){
        return {
            name: "NewField",
            description: "",
            type: $scope.PRIMITIVES[0],
            example: "",
            constraints: [],
            access: "PUBLIC"
        };
    };
    
    $scope.newConstraintDescription = function (){
        return {
            type: "regex",
            value: ""
        };
    };

    $scope.lookupPrimitiveJavaType = {
        "string" : "java.lang.String",
        "number" : "java.lang.Long",
        "boolean" : "java.lang.Boolean",
        "object" : "java.util.HashMap"
    };
        
    

    $scope.submit = function (schemaDescription){
        
        schemaDescription.fields.map(function(field){
            if ($scope.PRIMITIVES.indexOf(field.type) != -1) {
                field.reference = false;
                field.javaType = $scope.lookupPrimitiveJavaType[field.type.name];
                field.type = field.type.name;
            }
            else {
                field.reference = true;
                field.javaType = field.schema.javaType;
            }
            if (field.constraints.length > 0) {
                field.constraints = {pattern: {regex: field.constraints[0].value}};
            } else {
                field.constraints = null;
            }
        });


        schemaDescription.schema = "ClothoSchema";
        schemaDescription.language = "JSONSCHEMA";
        
        if (schemaDescription.superClass) schemaDescription.superClass = schemaDescription.parentSchema.id;
        //XXX: handle better
        delete schemaDescription.parentSchema;

        Clotho.create(schemaDescription).then(function (data){
            Clotho.get(data).then(function (data){
                $scope.schemas.push(data);
                $scope.created = $scope.newDescription();
            })
        });
    };

    $scope.create = function (schema) {
        var skeleton = {"schema":schema.name};
        Clotho.create(skeleton).then(function (id){
            $location.path("#/editor");
        });
    };

    $scope.getTypes = function (){
        return $scope.PRIMITIVES.concat($scope.schemas);
    };

    $scope.ACCESS_LIST = ["PRIVATE", "PUBLIC", "READONLY"];
    $scope.PRIMITIVES = [{name:"string"}, {name:"number"}, {name:"boolean"}, {name:"object"}];
    $scope.CONSTRAINTS = ["regexp"];

        Clotho.query({"schema":"Schema"}).then(function(data){
            $scope.schemas = data;
        });
    
    $scope.setCurrent($scope.nav[0]);
    $scope.created = $scope.newDescription();
}]);

