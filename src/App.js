import Home from './pages/Home/Home';
import Input from './pages/Input/Input';
import Prediction from './pages/Prediction/Prediction';
import {
  BrowserRouter as Router,
  Routes,
  Route,
  Link
} from 'react-router-dom';
import './App.css';

function App() {
  return (
    <Router>
      <Routes>
        <Route exact path='/' element={ <Home /> } />
        <Route exact path='/input' element={ <Input /> } />
        <Route exact path='/prediction' element={ <Prediction /> } />
      </Routes>
    </Router>
  );
}

export default App;
