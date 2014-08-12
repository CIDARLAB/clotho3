/**
 * @ngdog directive
 * @module clotho.commandbar
 *
 * @name clothoFunctionExecutor
 *
 * @attr function {Object} Function to use for executor
 * @attr functionResult {*=} Bind to result of run
 * @attr onExecute {Function=} Run when function is executed, with param $result
 */

angular.module('clotho.commandbar')
.directive('clothoFunctionExecutor', function (Clotho, ClothoSchemas) {
	return {
		scope : {
			"function" : '=',
			"onExecute" : '&?'
		},
		templateUrl: 'views/_command/executor.html',
		link: function (scope, element, attrs) {

			scope.isPrimitiveField = ClothoSchemas.isPrimitiveField;
			scope.schemaReadable = ClothoSchemas.mapSchemaIdToName;

			//future - filter
			scope.capitalize = function (word) {
				return word.substring(0, 1).toUpperCase() + word.substr(1);
			};

			function resetFunctionArgs() {
				scope.functionArgs = {};
			}

			function flattenArgs(functionObject, argsObject) {
				var arr = [];
				angular.forEach(functionObject.args, function (arg) {
					arr.push(argsObject[arg.name]);
				});
				return arr;
			}

			scope.executeFunction = function () {
				Clotho.run(scope.function.id, flattenArgs(scope.function, scope.functionArgs))
				.then(function (result) {
					scope.functionResult = result;
					scope.onExecute({$result : result});
				});
			};

			scope.clearArguments = function () {
				resetFunctionArgs();
				scope.functionResult = null;
			};

			scope.setArgument = function (name, id) {
				Clotho.get(id, {mute: true}).then(function (r) {
					scope.functionArgs[name] = r;
				});
			};

			scope.$watch('function.id', function (newid) {
				resetFunctionArgs();
			});
		}
	}
});