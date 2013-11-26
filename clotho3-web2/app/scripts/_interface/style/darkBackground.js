angular.module('clotho.interface')
  .directive('darkBackground', function ($document) {
    return {
      restrict: 'A',
      link: function postLink(scope, element, attrs) {
        //todo - convert to class, use toggleClass(), remove on scope.$destroy() or $location change
	      var oldBackground = $document.find('body').css('background');
				$document.find('body').css({background: '#eeeeee'});

	      scope.$on("$destroy", function() {
		      $document.find('body').css({background: oldBackground});
	      })
      }
    };
  });
