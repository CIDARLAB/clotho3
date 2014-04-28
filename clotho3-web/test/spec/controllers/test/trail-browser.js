'use strict';

describe('Controller: TestTrailBrowserCtrl', function () {

  // load the controller's module
  beforeEach(module('clothoApp'));

  var TestTrailBrowserCtrl,
    scope;

  // Initialize the controller and a mock scope
  beforeEach(inject(function ($controller, $rootScope) {
    scope = $rootScope.$new();
    TestTrailBrowserCtrl = $controller('TestTrailBrowserCtrl', {
      $scope: scope
    });
  }));

  it('should attach a list of awesomeThings to the scope', function () {
    expect(scope.awesomeThings.length).toBe(3);
  });
});
