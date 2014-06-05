angular.module('clotho.construction')
/**
 * @ngdoc directive
 * @name constructionFileView
 *
 * @restrict A
 *
 * @description
 * Visualizer for construction files.
 */
	.directive('constructionFile', function (ConstructionSimulator, clothoTokenCollectionFactory) {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionFile.html',
			scope : {
				inputFile : '=constructionFile',
				noProcess : '='
			},
			link : function (scope, element, attrs) {
				scope.$watch('inputFile', function (newfile) {
					scope.file = newfile;
				});

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