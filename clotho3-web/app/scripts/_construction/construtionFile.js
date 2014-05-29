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
	.directive('constructionFile', function (ConstructionSimulator) {
		return {
			restrict : 'A',
			templateUrl : 'views/_construction/constructionFile.html',
			scope : {
				file : '=constructionFile'
			},
			link : function (scope, element, attrs) {
				scope.$watch('file', function (newfile) {
					ConstructionSimulator.process(newfile).then(function (fileResult) {
						scope.computed = fileResult;
					});
				});
			}
		}
	});