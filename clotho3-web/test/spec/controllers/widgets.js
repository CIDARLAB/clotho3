'use strict';

describe('Controller: WidgetsCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var WidgetsCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    WidgetsCtrl = $controller('WidgetsCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
