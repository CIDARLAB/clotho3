'use strict';

Application.Directives.directive('sharableEditorOld', ['Collector', function(Collector) {
    /*
    //FUTURE - PROBLEMS
    - RESET --- form dirty not set when each input is a directive - need to emit or something
    - missing fields should be accounted for now
     */

    //TODO - load before page load, or have promise resolution re-trigger form generation
    //NOTE - complicated by directive having isolate scope -- need to inform
    var schemas = {};
    //schemas.Institution = Collector.retrieveModel('schema_institution');
    //schemas.Person = Collector.retrieveModel('schema_person');


    // FIXME -- multiple page collector doesn't work - see below (b/c using promises) --- need consistency in passing around models

    return {
        replace: false,
        require: 'form',
        restrict: 'A',
        controller: function($scope, $element, $attrs, Collector) {
            //sharable is inherited from parent, though could define it explicitly...
            $scope.schemaName = $attrs.schema;

            Collector.retrieveModel($scope.schemaName).then(function(result) {
                schemas[$scope.schemaName] = result;
                $scope.schema = result.schema;
                $scope.schema_custom = result.custom;

                //testing
                //console.log("schema promise fulfilled");
                //console.log($scope.schema);
            });

            //note -- need to access form constructor like: $scope[$scope.formName].$setPristine()
            $scope.formName = $attrs.name;

            $scope.editMode = false;
            $scope.formDirty = false;

            $scope.$watch('uuid', function(newval, oldval) {
                // note - won't change if explicity set in directive (due to prototypical inheritance)
                // testing:
                // console.log("SHAR_ED_DIR\tuuid changed: " + newval);
            });


            /*Collector.watch("inst_first", function(newData) {
                console.log("SHAR_ED_DIR\tCollector.watch: " + $scope.uuid);
                $scope.sharable = newData
            });*/

            //switch to 'edit' mode
            $scope.edit = function() {
                $scope.editMode = true;
            };

            //discard edits
            $scope.reset = function() {
                //FUTURE - if missing fields and invalid (e.g. email) won't get reset -- replace whole sharable object... or, check for $dirty fields
                Collector.retrieveModel($scope.uuid).then(function(result) {
                    $scope.sharable = result;
                });
                $scope[$scope.formName].$setPristine();
            };

            //save edits, switch to 'view'
            $scope.save = function() {
                //FUTURE - make sure we want this check -- don't want when elements have own scope
                /*
                if ($scope.formDirty) {
                    Collector.storeModel($scope.uuid, $scope.sharable);
                }
                */
                console.log("SHAR_ED_DIR\tstoring $scope.sharable");
                console.log($scope.sharable);
                Collector.storeModel($scope.uuid, $scope.sharable);
                $scope.editMode = false;
            };

            //discard edits, switch to 'view'
            $scope.discard = function() {
                $scope.reset();
                $scope.editMode = false;
            };

            //testing
            $scope.logScope = function() {
                console.log($scope);
            };

            //can't pass in the object (in current Angular version) or breaks before promise fulfilled, so:
            // true at the end compares values, omitting (false, default) compares references...
            // ...(and won't fire for watching object if only a child field changes)
            $scope.$watch('sharable', function(newVal, oldVal) {
                //console.log($scope.sharable);
            }, true);
        },
        compile: function compile(tElement, tAttrs, transclude) {

            return {
                pre: function preLink(scope, iElement, iAttrs, controller) {

                },
                post: function postLink(scope, iElement, iAttrs, controller) {

                    //todo - move to native angular code to handle this
                    //avoid complex statements in watch statements (e.g. do we need the else if?)
                    //this statement only works for when form elements do not have new scope
                    scope.$watch(iAttrs.name + '.$dirty', function (newValue, oldValue) {
                        if (newValue != oldValue && newValue === true) {
                            scope.formDirty = true;
                        } else if (newValue != oldValue && newValue === false) {
                            scope.formDirty = false;
                        } else {}
                    });
                }
            }
        }
    }
}]);
