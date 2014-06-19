angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionFileView
 *
 * @restrict A
 *
 * @description
 * Visualizer / Simulator GUI for construction files. Automatically preprocess by default (interpolate inputs on clotho)
 *
 * @attr constructionFile File to simulate /display
 * @attr noProcess Hide process button
 * @attr noPreprocess Do not automatically preprocess
 * @attr auto-process Automatically process
 *
 * @example
 *
 * <div construction-file="myConstructionFile" cf-preprocess no-process="true"></div>
 */
	.directive('constructionFile', function (ConstructionSimulator) {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionFile.html',
			scope : {
				inputFile : '=constructionFile',
				noProcess : '='
			},
			link : function (scope, element, attrs) {

				var autoPreprocess = angular.isUndefined(attrs.noPreprocess);
				var autoProcess = angular.isDefined(attrs.autoProcess);

				function resetProcessFlags () {
					angular.forEach(['processing', 'processed', 'processError'], function (flag) {
						scope[flag] = false;
					});
				}

				scope.$watch('inputFile', function (newfile) {
					resetProcessFlags();
					scope.file = newfile;
					if (autoProcess) {
						scope.process();
					}
					else if (autoPreprocess) {
						scope.preprocess();
					}
				});

				//todo - check all inputs
				scope.preprocess = function () {
					ConstructionSimulator.interpolateInputs(scope.file);
				};

				scope.process = function () {
					scope.processing = true;
					ConstructionSimulator.process(scope.file)
					.then(function onCFSuccess(processed) {
						scope.processed = true;
					}, function onCFError(partialFile) {
						scope.processError = true;
					})
					.finally(function () {
						scope.processing = false;
					});
				}
			}
		}
	});